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
//#define SMALL_TIMER

//#define PHASE_TIMER

size_t setup_grain_size = 20;
size_t update_degree_grain_size = 15;

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
inline void writeAdd(ET *a, ET b, size_t k) {
  volatile ET newV, oldV;
  do {oldV = *a; newV = oldV + b;}
  while (oldV > k && !CAS(a, oldV, newV));
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

#ifdef SMALL_TIMER
  Timer small_timer;
  size_t threshold = 20;
  double total_small_time = 0;
#endif

#ifdef PHASE_TIMER
  Timer initial_bin_setup_timer;
  Timer degree_update_phase_timer;
  Timer find_smallest_bin_phase_timer;
  double total_degree_update_time;
  double total_find_smallest_bin_time;
#endif

  total_t.Start();
  
  //t.Start();
 //start from the first round, and proceed

  // two element arrays for double buffering curr=iter&1, next=(iter+1)&1
  size_t shared_indexes[2] = {0, kMaxBin};
  size_t frontier_tails[2] = {0, 0};
  size_t first_frontier_tail = 0;
  //set up a boolean vector to keep track wether each vertex has been processed already

#ifdef PHASE_TIMER 
  //#pragma omp single 
    {
      initial_bin_setup_timer.Start();
    }
#endif

  pvector<bool> processed(g.num_nodes(), false);

  size_t start_bin_index = kMaxBin;
  

  #pragma omp parallel 
  {

    //doing a first pass to put every node into the right initial bin
    vector<vector<NodeID>> local_bins(1);



    #pragma omp for nowait schedule(dynamic, setup_grain_size)
    for (NodeID i = 0; i < g.num_nodes(); i++){
      size_t dest_bin = degree[i];
      //if (dest_bin < min_degree_threshold){
        if (dest_bin >= local_bins.size()){
          local_bins.resize(dest_bin+1);
        }//end of resize if
#ifdef DEBUG_ATOMICS
	  cout << "  putting node: " << i << " with degree: " << degree[i] <<  " into bin: " << dest_bin << endl;
#endif

	local_bins[dest_bin].push_back(i);

	//} // end of if less than min_degree_threshold
  }//end of for loop
    

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
    

 #ifdef PHASE_TIMER
    #pragma omp single 
    {
      initial_bin_setup_timer.Stop();
    }
 #endif

    size_t iter = 0;

    
    while (shared_indexes[iter & 1] != kMaxBin){
      size_t &curr_bin_index = shared_indexes[iter&1];
      size_t &next_bin_index = shared_indexes[(iter+1)&1];
      size_t &curr_frontier_tail = frontier_tails[iter&1];
      size_t &next_frontier_tail = frontier_tails[(iter+1)&1];
      size_t k = curr_bin_index;
 
#ifdef PROFILE     
      #pragma omp single
      {
        cout << "iter: " << iter << endl;
	cout << "k: " << k << endl;
        cout << " current frontier size:  " << curr_frontier_tail << endl;
	cout << " current bin index: " << curr_bin_index << endl;
      }
#endif
      
     /** if (curr_frontier_tail < 100){
	#pragma omp single 
	{

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
		size_t latest_degree = degree[ngh] - 1;
		//insert into the right bucket
		size_t dest_bin = latest_degree;
		if (dest_bin >= local_bins.size()){
		  local_bins.resize(dest_bin+1);
		}
		local_bins[dest_bin].push_back(ngh);
	      }//end of if degree[ngh] > k
	    }//end of inner for
	  }//end of outer for
	}//end of pragma omp si 
	}//end of the serial version 
	else {**/

#ifdef SMALL_TIMER
  #pragma omp single 
  {
    if ( curr_frontier_tail < threshold){small_timer.Start();}
  }
#endif

#ifdef PHASE_TIMER
  #pragma omp single 
  {
     degree_update_phase_timer.Start();
  }
#endif

      #pragma omp for schedule (dynamic, update_degree_grain_size)
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
	    //NodeID latest_degree = fetch_and_add(degree[ngh],-1) - 1;
    writeAdd(&degree[ngh],-1, k);
	    NodeID latest_degree = degree[ngh];

#ifdef DEBUG_ATOMICS
	    cout << "  node: " << ngh << " with latest degree: " << latest_degree << endl;
#endif

	    //if (latest_degree >= k){ //only update if it is more than the k degree
	      //insert into the right bucket
	    NodeID dest_bin = latest_degree;
	      if (dest_bin >= local_bins.size()){
		local_bins.resize(dest_bin+1);
	      }
#ifdef DEBUG_ATOMICS
	      cout << "  push node: " << ngh << " into bin: " << dest_bin << endl;
#endif
	      local_bins[dest_bin].push_back(ngh);
	    } 
	  //}//end of if statement
        } //end of inner for
      }//end of outer for

#ifdef PHASE_TIMER
      #pragma omp single 
      {
        degree_update_phase_timer.Stop();
	total_degree_update_time += degree_update_phase_timer.Seconds();
	find_smallest_bin_phase_timer.Start();
      }
#endif

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


#ifdef SMALL_TIMER
     #pragma omp single 
     {
      if (curr_frontier_tail < threshold){
	small_timer.Stop();
        total_small_time += small_timer.Millisecs();
	//cout << "small timer seconds: " << small_timer.Millisecs() << endl;
      }
     }
#endif


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

#ifdef PHASE_TIMER
     #pragma omp single 
     {
        find_smallest_bin_phase_timer.Stop();
	total_find_smallest_bin_time += find_smallest_bin_phase_timer.Seconds();
      }
#endif

    }//end of while
    #pragma omp single
    cout << "number of iter: " << iter << endl;
  }//end of the parallel region
  
  total_t.Stop();
  cout << "total exec time: " << total_t.Seconds() << endl;

#ifdef SMALL_TIMER
  cout << "small exec time: " << total_small_time/1000 << endl;
#endif

#ifdef PHASE_TIMER
cout << " initial bin setup time: " << initial_bin_setup_timer.Millisecs()/1000 << endl;
cout << " degree update time: " << total_degree_update_time  << endl;
cout << " find next bin time: " << total_find_smallest_bin_time   <<endl;
#endif

  int max_core = 0;
  int sum_core = 0;
  for (NodeID i = 0; i < g.num_nodes(); i++){
    if (degree[i] > max_core) max_core = degree[i];
    sum_core += degree[i];
#ifdef PRINT_CORES
    cout << "Node " <<  i << " core: " <<  degree[i] << endl;
#endif
  }

  cout << "max of core: " << max_core << endl;
  cout << "sum of core: " << sum_core << endl;
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
