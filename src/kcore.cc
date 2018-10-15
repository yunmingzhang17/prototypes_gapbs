// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <algorithm>
#include <iostream>
#include <vector>

#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "pvector.h"

#define DEBUG_ATOMICS
#define TIME_ATOMICS

using namespace std;

const size_t kMaxBin = numeric_limits<size_t>::max()/2;

pvector<NodeID> kcore_atomics (const Graph &g){


  pvector<NodeID> degree(g.num_nodes());
  std::cout << "running kcore" << std::endl;
  //iterate though the current bucket 
#ifdef DEBUG_ATOMICS
  size_t deg_1_count = 0;
#endif

  for (NodeID n : g.vertices()){
    degree[n] = g.out_degree(n);
    #ifdef DEBUG_ATOMICS
    if (degree[n] == 1) deg_1_count++;
    #endif
  }

  
#ifdef DEBUG_ATOMICS
  cout << "actual deg 1 vertices count: " << deg_1_count <<  endl;
#endif

  pvector<NodeID> frontier(g.num_edges_directed());


  Timer t;
  Timer total_t;

  total_t.Start();
  
  t.Start();

  //do something to initilize the frontier with degree 1 nodes efficiently 
  size_t frontier_tail = 0;

  #pragma omp parallel 
  {
    vector<NodeID> local_bin(0);
    #pragma omp for nowait schedule(dynamic, 64)
    for (NodeID i = 0; i < g.num_nodes(); i++){
      if (degree[i] == 1){
	local_bin.push_back(i);
      }
    }
    #pragma omp barrier 

    size_t copy_start = fetch_and_add(frontier_tail, 
				      local_bin.size());
    copy(local_bin.begin(), local_bin.end(), frontier.data() + copy_start);
    local_bin.resize(0);
    #pragma omp barrier
    
  }

 t.Stop();

#ifdef TIME_ATOMICS
 cout << "first round took: " << t.Millisecs()/1000 << endl;
#endif

#ifdef DEBUG_ATOMICS
  std::cout << "measured number of degree 1 vertices: " << frontier_tail << endl;
#endif  
  
  //start from the first round, and proceed

  // two element arrays for double buffering curr=iter&1, next=(iter+1)&1
  size_t shared_indexes[2] = {0, kMaxBin};
  size_t frontier_tails[2] = {frontier_tail, 0};
  
  //set up a boolean vector to keep track wether each vertex has been processed already
  pvector<bool> processed(g.num_nodes(), false);
  
t.Start();



  #pragma omp parallel
  {
    size_t iter = 0;
    vector<vector<NodeID>> local_bins(0);
    
    while (shared_indexes[iter & 1] != kMaxBin){
      size_t &curr_bin_index = shared_indexes[iter&1];
      size_t &next_bin_index = shared_indexes[(iter+1)&1];
      size_t &curr_frontier_tail = frontier_tails[iter&1];
      size_t &next_frontier_tail = frontier_tails[(iter+1)&1];
      size_t k = curr_bin_index + 2;

      #pragma omp single
      cout << " current frontier size:  " << curr_frontier_tail << endl;

      #pragma omp for nowait schedule (dynamic, 64)
      for (size_t i = 0; i < curr_frontier_tail; i++){
        NodeID u = frontier[i];
	// if the node is already processed in an earlier bin
	// then skip the current node
	if (processed[u]) continue;
	// set the node to be processed after this round (removed)
	else processed[u] = true;

	for (NodeID ngh : g.out_neigh(u)){
          if (degree[ngh] > k) {
	    //update the degree of the neighbor
	    // this value should be unique across threads, no duplicated vertices in the same bin
	    size_t latest_degree = fetch_and_add(degree[ngh],-1);

	    if (latest_degree >= k){ //only update if it is more than the k degree
	      //insert into the right bucket
	      size_t dest_bin = latest_degree;
	      if (dest_bin >= local_bins.size()){
		local_bins.resize(dest_bin+1);
	      }
	      local_bins[dest_bin].push_back(ngh);
	    } 
	  }//end of if statement
        } //end of inner for
      }//end of outer for

      for (size_t i = curr_bin_index; i < local_bins.size(); i++){
	if (!local_bins[i].empty()){
	  #pragma omp critical
	  next_bin_index = min(next_bin_index, i);
	  break;
	}
      }// end of for loop to find the next bin index



      #pragma omp barrier
      
      //cout << "next_bin_index:" << next_bin_index << endl;

      #pragma omp single nowait
      {
      //t.Stop();
      // PrintStep(curr_bin_index, t.Millisecs(), curr_frontier_tail);
      // t.Start();
        curr_bin_index = kMaxBin;
        curr_frontier_tail = 0;
      }

      if (next_bin_index < local_bins.size()) {
        size_t copy_start = fetch_and_add(next_frontier_tail,
                                          local_bins[next_bin_index].size());
        copy(local_bins[next_bin_index].begin(),
             local_bins[next_bin_index].end(), frontier.data() + copy_start);
        local_bins[next_bin_index].resize(0);
      }
      iter++;
      #pragma omp barrier

    }//end of while
    #pragma omp single
    cout << "number of iter: " << iter << endl;
  }//end of the parallel region
  
  total_t.Stop();
  cout << "total exec time: " << total_t.Millisecs()/1000 << endl;

  int max_core = 0;
#pragma omp parallel for reduction(max: max_core)
  for (size_t i = 0; i < g.num_nodes(); i++){
    max_core = degree[i];
  }

  cout << "max of core: " << max_core << endl;

  return degree;
}

int main(int argc, char* argv[]) {
  CLApp cli(argc, argv, "kcore");
  if (!cli.ParseArgs())
    return -1;
  Builder b(cli);
  Graph g = b.MakeGraph();
  std::cout << "num vertices: " << g.num_nodes() << std::endl;
  for (int trail = 0; trail < 3; trail++){
    kcore_atomics(g);
  }
  return 0;
}
