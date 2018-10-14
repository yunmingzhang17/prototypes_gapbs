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
  // two element arrays for double buffering curr=iter&1, next=(iter+1)&1
  size_t shared_indexes[2] = {0, kMaxBin};
  size_t frontier_tails[2] = {1, 0};

  Timer t;
  
  t.Start();

  //do something to initilize the frontier with degree 0 nodes efficiently 
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

#ifdef DEBUG_ATOMICS
  std::cout << "measured number of degree 1 vertices: " << frontier_tail << endl;
#endif  
  
  //start from the first round, and proceed

  
  t.Stop();

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
