// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <algorithm>
#include <iostream>
#include <vector>
#include <cinttypes>
#include <queue>
#include <limits>

#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "pvector.h"

#include <unordered_map>


//#define DEBUG_ATOMICS
#define TIME_ATOMICS
//#define PRINT_CORES
//#define PROFILE



using namespace std;




template <class ET>
inline bool CAS(ET *ptr, ET oldv, ET newv) {
  if (sizeof(ET) == 1) {
    return __sync_bool_compare_and_swap((bool*)ptr, *((bool*)&oldv), *((bool*)&newv));
  } else if (sizeof(ET) == 4) {
    return __sync_bool_compare_and_swap((int*)ptr, *((int*)&oldv), *((int*)&newv));
  } else if (sizeof(ET) == 8) {
    return __sync_bool_compare_and_swap((long*)ptr, *((long*)&oldv), *((long*)&newv));
  }
  else {
    std::cout << "CAS bad length : " << sizeof(ET) << std::endl;
    abort();
  }
}


template <class ET>
inline void writeAdd(ET *a, ET b) {
  volatile ET newV, oldV;
  do {oldV = *a; newV = oldV + b;}
  while (!CAS(a, oldV, newV));
}

const size_t kMaxBin = numeric_limits<size_t>::max()/2;

bool has_unprocessed(vector<NodeID> bin, pvector<bool> &processed){
  for (NodeID i : bin){
    if (processed[i] == false) return true;
  }
  return false;
}

NodeID* kcore_atomics (const Graph &g){


  NodeID* degree = new NodeID[g.num_nodes()];
  std::cout << "running kcore" << std::endl;
  //iterate though the current bucket 
#ifdef DEBUG_ATOMICS
  //size_t deg_1_count = 0;
#endif

  //size_t min_degree = numeric_limits<size_t>::max()/2;

#pragma omp parallel for
  for (NodeID i = 0; i < g.num_nodes(); i++){
    degree[i] = g.out_degree(i);
    //if (min_degree > degree[n] && degree[n] != 0) min_degree = degree[n];
  }

  
#ifdef DEBUG_ATOMICS
    //cout << "actual deg 1 vertices count: " << deg_1_count <<  endl;
    //cout << "min degree: " << min_degree << endl;
#endif

  pvector<NodeID> frontier(g.num_edges_directed());


  Timer t;
  Timer total_t;

  total_t.Start();
  
  //t.Start();
 //start from the first round, and proceed

  // two element arrays for double buffering curr=iter&1, next=(iter+1)&1
  size_t shared_indexes[2] = {0, kMaxBin};
  size_t frontier_tails[2] = {0, 0};
  size_t first_frontier_tail = 0;
  //set up a boolean vector to keep track wether each vertex has been processed already
  pvector<bool> processed(g.num_nodes(), false);

  size_t start_bin_index = kMaxBin;


  #pragma omp parallel 
  {

    //doing a first pass to put every node into the right initial bin
    vector<vector<NodeID>> local_bins(1);


    unordered_map<NodeID, int> update_buffer;
    //vector<NodeID> update_keys;


    #pragma omp for nowait schedule(static, 64)
    for (NodeID i = 0; i < g.num_nodes(); i++){
      size_t dest_bin = degree[i];
      if (dest_bin >= local_bins.size()){
        local_bins.resize(dest_bin+1);
      }

#ifdef DEBUG_ATOMICS
      cout << "  putting node: " << i << " with degree: " << degree[i] <<  " into bin: " << dest_bin << endl;
#endif

      local_bins[dest_bin].push_back(i);

    }

    //find the bin to start (the smallest degree)

    for (size_t i = 0; i < local_bins.size(); i++){
      if (!local_bins[i].empty()){
        #pragma omp critical
	start_bin_index = min(start_bin_index, i);
	break;
      }
    }// end of for loop to find the smallest bin index

    #pragma omp barrier 
    
    size_t copy_start = fetch_and_add(first_frontier_tail, 
				      local_bins[start_bin_index].size());
    copy(local_bins[start_bin_index].begin(), local_bins[start_bin_index].end(), frontier.data() + copy_start);
    //release the  bin after being copied
    local_bins[start_bin_index].resize(0);    

#ifdef DEBUG_ATOMICS
    cout << "start bin index: " << start_bin_index << endl;
    cout << "start bin size: " << first_frontier_tail << endl;
#endif

    #pragma omp single
    {
      shared_indexes[0] = start_bin_index;
      frontier_tails[0] = first_frontier_tail;
    }
    #pragma omp barrier
    

 

    size_t iter = 0;

    
    while (shared_indexes[iter & 1] != kMaxBin){
      size_t &curr_bin_index = shared_indexes[iter&1];
      size_t &next_bin_index = shared_indexes[(iter+1)&1];
      size_t &curr_frontier_tail = frontier_tails[iter&1];
      size_t &next_frontier_tail = frontier_tails[(iter+1)&1];
      size_t k = curr_bin_index;


      update_buffer.clear();
      //update_keys.resize(0);
 
#ifdef PROFILE     
      #pragma omp single
      {
        cout << "iter: " << iter << endl;
	cout << "k: " << k << endl;
        cout << " current frontier size:  " << curr_frontier_tail << endl;
	cout << " current bin index: " << curr_bin_index << endl;
      }
#endif
      

      #pragma omp for nowait schedule (dynamic, 10)
      for (size_t i = 0; i < curr_frontier_tail; i++){
        NodeID u = frontier[i];
	// if the node is already processed in an earlier bin
	// then skip the current node
	if (processed[u]) continue;
	// set the node to be processed after this round (removed)
	else processed[u] = true;

#ifdef DEBUG_ATOMICS
	cout << "current node: " << u << endl;
#endif

	for (NodeID ngh : g.out_neigh(u)){

#ifdef DEBUG_ATOMICS
          cout << "  ngh " << ngh <<  " with degree: " << degree[ngh] << endl;
#endif

          if (degree[ngh] > k) {
	    //update the degree of the neighbor
	    // this value should be unique across threads, no duplicated vertices in the same bin
#ifdef DEBUG_ATOMICS
	    cout << "reducing local buffered delta by -1 for node: " << ngh << endl;
#endif
	    if (update_buffer.find(ngh) != update_buffer.end()){
	      update_buffer[ngh] += -1;
	    } else {
	      update_buffer[ngh] = -1;
	      //update_keys.push_back(ngh);
	    }
	  }//end of if statement
	} //end of inner for
      }//end of outer for


      for (unordered_map<NodeID, int>:: iterator itr = update_buffer.begin(); itr != update_buffer.end(); itr++) {
      //    #pragma omp parallel for
      //for (int i = 0; i < update_keys.size(); i++){ 
	
	NodeID node = itr->first;
	int delta = itr->second;
	writeAdd(&degree[node],delta);
	NodeID latest_degree = degree[node] > curr_bin_index ? degree[node]: curr_bin_index;
	NodeID dest_bin = latest_degree;
#ifdef DEBUG_ATOMICS
	cout << "  node: " << node << " with latest degree: " << latest_degree << " delta: " << delta  << endl;
#endif
	if (dest_bin >= local_bins.size()){
	  local_bins.resize(dest_bin+1);
	}
#ifdef DEBUG_ATOMICS
	cout << "  push node: " << node << " into bin: " << dest_bin << endl;
#endif
	local_bins[dest_bin].push_back(node);
      } 

	


      for (size_t i = curr_bin_index; i < local_bins.size(); i++){
      //cout << "index: " << i << endl;
      //	cout << "next bin index: " << next_bin_index << endl;
	if (!local_bins[i].empty() && has_unprocessed(local_bins[i], processed)){
	  #pragma omp critical
	  {

    //#ifdef PROFILE
	   
	   //#endif
	   next_bin_index = min(next_bin_index, i);
	  }
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
  for (NodeID i = 0; i < g.num_nodes(); i++){
    if (degree[i] > max_core) max_core = degree[i];
#ifdef PRINT_CORES
    cout << "Node " <<  i << " core: " <<  degree[i] << endl;
#endif
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
  std::cout << "num edges: " << g.num_edges() << std::endl;
  for (int trail = 0; trail < 5; trail++){
    kcore_atomics(g);
  }
  return 0;
}
