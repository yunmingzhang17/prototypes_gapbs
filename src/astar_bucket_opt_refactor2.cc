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

Kernel: Point to Point Shortest Path (PPSP)
Author: Yunming Zhang (experimental)

Returns array of distances for all vertices from given source vertex

*/


using namespace std;


extern void load_coords(std::string filename, int num_nodes);
extern double calculate_distance (NodeID source, NodeID destination);

const WeightT kDistInf = numeric_limits<WeightT>::max()/2;
//const size_t kMaxBin = numeric_limits<size_t>::max()/2;
const size_t bin_size_threshold = 500;
EagerPriorityQueue<WeightT> * pq;
WeightT* __restrict dist_array;
WeightT* __restrict est_dist_array;
WeightT input_delta;
NodeID dest;

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
inline bool writeMin(ET *a, ET b) {
  ET c; bool r=0;
  do c = *a; 
  while (c > b && !(r=CAS(a,c,b)));
  return r;
}

struct edge_update_func
{
  void operator()(vector<vector<NodeID> >& local_bins, NodeID src, NodeID dst, WeightT wt){
    WeightT new_dist = dist_array[src] + wt;
    //bool changed = writeMin(&dist_array[dst], new_dist);

    if (new_dist < dist_array[dst]) dist_array[dst] = new_dist;

    //if (changed){
      WeightT old_est_dist = est_dist_array[dst];
      WeightT new1 = dist_array[dst] + calculate_distance(dst, dest);
      //cout << "src node: " << src << " dst: " << dst << " wt: " << wt  << endl;
      //cout << "dst node: " << dst <<" new est dist: " << new1 << endl;
      //cout << "dst node: " << dst <<" actual dist: " << dist_array[dst] << endl;
      WeightT new2 = est_dist_array[dst];
      WeightT new_est_dist = new1 > new2 ? new1 : new2;
      if (new_est_dist < old_est_dist)
	update_priority_min<WeightT>()(pq, local_bins, dst, old_est_dist, new_est_dist);
      //}
  }
};

struct src_filter_func
{
  bool operator()(NodeID v){
    return (est_dist_array[v] < est_dist_array[dest] && est_dist_array[v] >= input_delta * static_cast<WeightT>(pq->get_current_priority()));
    //return est_dist_array[v] >= input_delta * static_cast<WeightT>(pq->get_current_priority());
  }
};

struct while_cond_func
{
  bool operator()(){
    cout << "estimated distance to dest: " << est_dist_array[dest] << endl;
    cout << "current priority: " << pq->get_current_priority() << endl;

    //return est_dist_array[dest]/input_delta >= pq->get_current_priority();
    return !pq->finished();
  }
};


pvector<WeightT> PPDeltaStep(const WGraph &g, NodeID source, NodeID dest, WeightT delta) {
  Timer t;
  pvector<WeightT> dist(g.num_nodes(), kDistInf);
  t.Start();

  //resetng the distances                                                                             
  #pragma omp parallel for
  for (int i = 0; i < g.num_nodes(); i++){
    dist_array[i] = kDistInf;
    est_dist_array[i] = kDistInf;
  }

  dist_array[source] = 0;
  est_dist_array[source] = calculate_distance(source, dest);
  
  OrderedProcessingOperatorWithMerge(pq, g, dist_array, src_filter_func(), while_cond_func(), edge_update_func(), 1000,  source);

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
bool SSSPVerifier(const WGraph &g, NodeID source, NodeID dest, 
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
	//cout << "src: " << u << " dst: " << wn.v << " wt: " << wn.w << endl;
        if (td + wn.w < oracle_dist[wn.v]) {
          oracle_dist[wn.v] = td + wn.w;
          mq.push(make_pair(td + wn.w, wn.v));
        }
      }
    }
  }


  // Report any mismatches
  //bool all_ok = true;
  //for (NodeID n : g.vertices()) {
  //if (dist_to_test[n] != oracle_dist[n]) {
  //    cout << n << ": " << dist_to_test[n] << " != " << oracle_dist[n] << endl;
  //    all_ok = false;
  //  }
  //}

  bool all_ok = false;
  if (dist_to_test[dest] == oracle_dist[dest]) all_ok = true;
  else {
    cout << "measured dist: " << dist_to_test[dest] << endl;
    cout << "oracle dist: " << oracle_dist[dest] << endl;
  }
  
  return all_ok;
}



int main(int argc, char* argv[]) {
  CLPPSP<WeightT> cli(argc, argv, "AStar search");
  if (!cli.ParseArgs())
    return -1;
  WeightedBuilder b(cli);
  WGraph g = b.MakeGraph();
  load_coords(cli.filename(), g.num_nodes());

  dist_array = new WeightT[g.num_nodes()];
  for (int i = 0; i < g.num_nodes(); i++){
    dist_array[i] = kDistInf;
  }

  est_dist_array = new WeightT[g.num_nodes()];
  for (int i = 0; i < g.num_nodes(); i++){
    est_dist_array[i] = kDistInf;
  }

  input_delta = cli.delta();
  pq = new EagerPriorityQueue<WeightT>(est_dist_array, input_delta);
  dest = cli.dst_vertex();
 
#ifdef ORIG

  SourcePicker<WGraph> sp(g, cli.start_vertex());
  auto SSSPBound = [&sp, &cli] (const WGraph &g) {
    // here we call the source picker two times, first for the source, then for the destination
    return PPDeltaStep(g, sp.PickNext(), sp.PickNext(),  cli.delta());
  };
  SourcePicker<WGraph> vsp(g, cli.start_vertex());
  auto VerifierBound = [&vsp] (const WGraph &g, const pvector<WeightT> &dist) {
    return SSSPVerifier(g, vsp.PickNext(), vsp.PickNext(), dist);
  };
  BenchmarkKernel(cli, g, SSSPBound, PrintSSSPStats, VerifierBound);
#else 

  if (cli.start_vertex() == -1) {
    cout << "invalide starting point -1" << endl;
    return 0;
  }

  auto SSSPBound = [ &cli] (const WGraph &g) {
    // here we call the source picker two times, first for the source, then for the destination
    return PPDeltaStep(g, cli.start_vertex(), cli.dst_vertex(),  cli.delta());
  };
  auto VerifierBound = [&cli] (const WGraph &g, const pvector<WeightT> &dist) {
    return SSSPVerifier(g, cli.start_vertex(), cli.dst_vertex(), dist);
  };

  int num_trails = 6;
  for (int i = 0; i < num_trails; i++)
    BenchmarkKernel(cli, g, SSSPBound, PrintSSSPStats, VerifierBound);


  cout << "starting point: " << cli.start_vertex() << endl;
  cout << "ending point: " << cli.dst_vertex() << endl;

#endif


  return 0;
}
