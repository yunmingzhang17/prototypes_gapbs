// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <cinttypes>
#include <limits>
#include <iostream>
#include <queue>
#include <vector>


#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "platform_atomics.h"
#include "pvector.h"
#include "timer.h"
#include "eager_priority_queue.h"
#include "ordered_processing.h"

/*
GAP Benchmark Suite
Kernel: Single-source Shortest Paths (SSSP)
Author: Scott Beamer

Returns array of distances for all vertices from given source vertex

This SSSP implementation makes use of the ∆-stepping algorithm [1]. The type
used for weights and distances (WeightT) is typedefined in benchmark.h. The
delta parameter (-d) should be set for each input graph.

The bins of width delta are actually all thread-local and of type std::vector
so they can grow but are otherwise capacity-proportional. Each iteration is
done in two phases separated by barriers. In the first phase, the current
shared bin is processed by all threads. As they find vertices whose distance
they are able to improve, they add them to their thread-local bins. During this
phase, each thread also votes on what the next bin should be (smallest
non-empty bin). In the next phase, each thread copies their selected
thread-local bin into the shared bin.

Once a vertex is added to a bin, it is not removed, even if its distance is
later updated and it now appears in a lower bin. We find ignoring vertices if
their current distance is less than the min distance for the bin to remove
enough redundant work that this is faster than removing the vertex from older
bins.

[1] Ulrich Meyer and Peter Sanders. "δ-stepping: a parallelizable shortest path
    algorithm." Journal of Algorithms, 49(1):114–152, 2003.
*/


using namespace std;

const WeightT kDistInf = numeric_limits<WeightT>::max()/2;


EagerPriorityQueue<WeightT> * pq; 
WeightT* __restrict dist_array;
WeightT input_delta;


struct edge_update_func
{
  void operator()(vector<vector<NodeID> >& local_bins, NodeID src, NodeID dst, WeightT wt){
  WeightT old_dist = dist_array[dst];
  WeightT new_dist = dist_array[src] + wt;
  update_priority_min<WeightT>()(pq, local_bins, dst, old_dist, new_dist);
  }
};

struct src_filter_func
{
  bool operator()(NodeID v){
    return (dist_array[v] >= input_delta * static_cast<WeightT>(pq->get_current_priority()));
  }
};

struct while_cond_func
{
  bool operator()(){
    return !pq->finished();
  }
};


pvector<WeightT> DeltaStep(const WGraph &g, NodeID source, WeightT delta) {
  Timer t;
  pvector<WeightT> dist(g.num_nodes(), kDistInf);

  t.Start();
  //resetng the distances
  #pragma omp parallel for
  for (int i = 0; i < g.num_nodes(); i++){
    dist_array[i] = kDistInf;
  }

  dist_array[source] = 0;
  
  //void OrderProcessingOperator(const WGraph &g, WeightT delta,Priority* dist_array, SrcFilter src_filter, WhileCond while_cond, NodeID optional_source_node){
  OrderedProcessingOperatorNoMerge(pq, g,  src_filter_func(), while_cond_func(), edge_update_func(),  source);

  t.Stop();
  cout << "DeltaStep took: " << t.Seconds() << endl;
  
  #pragma omp parallel for
  for (int i = 0; i < g.num_nodes(); i++){
    dist[i] = dist_array[i];
  }  

  return dist;
}


void PrintSSSPStats(const WGraph &g, const pvector<WeightT> &dist) {
  auto NotInf = [](WeightT d) { return d != kDistInf; };
  int64_t num_reached = count_if(dist.begin(), dist.end(), NotInf);
  cout << "SSSP Tree reaches " << num_reached << " nodes" << endl;
}


// Compares against simple serial implementation
bool SSSPVerifier(const WGraph &g, NodeID source,
                  const pvector<WeightT> &dist_to_test) {

  // Serial Dijkstra implementation to get oracle distances
  pvector<WeightT> oracle_dist(g.num_nodes(), kDistInf);
  oracle_dist[source] = 0;
  typedef pair<WeightT, NodeID> WN;
  priority_queue<WN, vector<WN>, greater<WN>> mq;
  mq.push(make_pair(0, source));
  while (!mq.empty()) {
    WeightT td = mq.top().first;
    NodeID u = mq.top().second;
    mq.pop();
    if (td == oracle_dist[u]) {
      for (WNode wn : g.out_neigh(u)) {
        if (td + wn.w < oracle_dist[wn.v]) {
          oracle_dist[wn.v] = td + wn.w;
          mq.push(make_pair(td + wn.w, wn.v));
        }
      }
    }
  }
  // Report any mismatches
  bool all_ok = true;
  for (NodeID n : g.vertices()) {
    if (dist_to_test[n] != oracle_dist[n]) {
      cout << n << ": " << dist_to_test[n] << " != " << oracle_dist[n] << endl;
      all_ok = false;
    }
  }
  return all_ok;
}


int main(int argc, char* argv[]) {
  CLDelta<WeightT> cli(argc, argv, "single-source shortest-path");
  if (!cli.ParseArgs())
    return -1;
  WeightedBuilder b(cli);
  WGraph g = b.MakeGraph();
  NodeID starting_node = cli.start_vertex();
  SourcePicker<WGraph> sp(g, starting_node);

  dist_array = new WeightT[g.num_nodes()];
  for (int i = 0; i < g.num_nodes(); i++){
    dist_array[i] = kDistInf;
  }
  
  input_delta = cli.delta();
  pq = new EagerPriorityQueue<WeightT>(dist_array, input_delta);
  

  auto SSSPBound = [&sp, &cli] (const WGraph &g) {
    return DeltaStep(g, sp.PickNext(), cli.delta());
  };


  SourcePicker<WGraph> vsp(g, cli.start_vertex());
  auto VerifierBound = [&vsp] (const WGraph &g, const pvector<WeightT> &dist) {
    return SSSPVerifier(g, vsp.PickNext(), dist);
  };
  //BenchmarkKernel(cli, g, SSSPBound, PrintSSSPStats, VerifierBound);

  int num_trails = 6;
  for (int i = 0; i < num_trails; i++)
    BenchmarkKernel(cli, g, SSSPBound, PrintSSSPStats, VerifierBound);

  cout << "starting node: " << starting_node << endl;
  cout << "starting node degree: " << g.out_degree(starting_node) << endl;
  return 0;
}
