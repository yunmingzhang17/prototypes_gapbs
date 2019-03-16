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

const WeightT kDistInf = numeric_limits<WeightT>::max()/2;
//const size_t kMaxBin = numeric_limits<size_t>::max()/2;
const size_t bin_size_threshold = 500;
EagerPriorityQueue<WeightT> * pq;
WeightT* __restrict dist_array;
WeightT input_delta;
NodeID dest;

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
    return (dist_array[v] < dist_array[dest] && dist_array[v] >= input_delta * static_cast<WeightT>(pq->get_current_priority()));
  }
};

struct while_cond_func
{
  bool operator()(){
    return dist_array[dest]/input_delta >= pq->get_current_priority();
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
  }

  dist_array[source] = 0;

  
  OrderedProcessingOperatorNoMerge(pq, g,  src_filter_func(), while_cond_func(), edge_update_func(), source);

  t.Stop();
  cout << "DeltaStep took: " << t.Seconds() << endl;

  #pragma omp parallel for
  for (int i = 0; i < g.num_nodes(); i++){
    dist[i] = dist_array[i];
  }

  return dist;

}


pvector<WeightT> PPDeltaStep_old(const WGraph &g, NodeID source, NodeID dest, WeightT delta) {
  Timer t;
  pvector<WeightT> dist(g.num_nodes(), kDistInf);
  dist[source] = 0;
  pvector<NodeID> frontier(g.num_edges_directed());
  // two element arrays for double buffering curr=iter&1, next=(iter+1)&1
  size_t shared_indexes[2] = {0, kMaxBin};
  size_t frontier_tails[2] = {1, 0};
  frontier[0] = source;
  //t.Start();
  #pragma omp parallel
  {
    vector<vector<NodeID> > local_bins(0);
    size_t iter = 0;
    while (shared_indexes[iter&1] != kMaxBin) {
      size_t &curr_bin_index = shared_indexes[iter&1];
      size_t &next_bin_index = shared_indexes[(iter+1)&1];
      size_t &curr_frontier_tail = frontier_tails[iter&1];
      size_t &next_frontier_tail = frontier_tails[(iter+1)&1];
      #pragma omp for nowait schedule(dynamic, 64)
      for (size_t i=0; i < curr_frontier_tail; i++) {
        NodeID u = frontier[i];

	//prune this point if it is not going to reach destination
	//it is possible that this point's distance would be reduced later
	//then it should be processed at a later time
	if (dist[u] > dist[dest]) continue;

        if (dist[u] >= delta * static_cast<WeightT>(curr_bin_index)) {
          for (WNode wn : g.out_neigh(u)) {
            WeightT old_dist = dist[wn.v];
            WeightT new_dist = dist[u] + wn.w;
            if (new_dist < old_dist) {
              bool changed_dist = true;
              while (!compare_and_swap(dist[wn.v], old_dist, new_dist)) {
                old_dist = dist[wn.v];
                if (old_dist <= new_dist) {
                  changed_dist = false;
                  break;
                }
              }
              if (changed_dist) {
                size_t dest_bin = new_dist/delta;
                if (dest_bin >= local_bins.size()) {
                  local_bins.resize(dest_bin+1);
                }
                local_bins[dest_bin].push_back(wn.v);
              }
            }
          }
        }
      }



      int round = 0;
      while (local_bins.size() > 0 && curr_bin_index < local_bins.size() && !local_bins[curr_bin_index].empty()){
      round++;
      size_t cur_bin_size = local_bins[curr_bin_index].size();
      if (cur_bin_size > bin_size_threshold) break;

      vector<NodeID> cur_bin_copy = local_bins[curr_bin_index];
      local_bins[curr_bin_index].resize(0);

      for (size_t i=0; i < cur_bin_size; i++) {
        NodeID u = cur_bin_copy[i];
        if (dist[u] >= delta * static_cast<WeightT>(curr_bin_index)) {
          for (WNode wn : g.out_neigh(u)) {
            WeightT old_dist = dist[wn.v];
            WeightT new_dist = dist[u] + wn.w;
            if (new_dist < old_dist) {
              bool changed_dist = true;
              while (!compare_and_swap(dist[wn.v], old_dist, new_dist)) {
                old_dist = dist[wn.v];
                if (old_dist <= new_dist) {
                  changed_dist = false;
                  break;
                }
              }
              if (changed_dist) {
                size_t dest_bin = new_dist/delta;
                if (dest_bin >= local_bins.size()) {
                  local_bins.resize(dest_bin+1);
                }
		local_bins[dest_bin].push_back(wn.v);

              }
            }
          }
        }
    }//end of outer for loop on curr_frontier_tail 

    }



      for (size_t i=curr_bin_index; i < local_bins.size(); i++) {
        if (!local_bins[i].empty()) {
          #pragma omp critical
          next_bin_index = min(next_bin_index, i);
          break;
        }
      }
      #pragma omp barrier
      
      //if next_bin_index is different from current_bin_index
      if (next_bin_index > curr_bin_index){
	// if the node is already processed
	if (dist[dest]/delta < next_bin_index) {
	  
	  #ifdef DEBUG_RESULT
	  cout << "src: " << source << endl;
	  cout << "dest: " << dest << endl;
	  cout << "current bin index: " << curr_bin_index << endl;
	  cout << "next bin index: " << next_bin_index << endl;
	  cout << "dest distance: " << dist[dest] << endl;
	  #endif

	  break;
	}
      }

      #pragma omp single nowait
      {
      //t.Stop();
      // PrintStep(curr_bin_index, t.Millisecs(), curr_frontier_tail);
      //  t.Start();
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
    }
    #pragma omp single
    cout << "took " << iter << " iterations" << endl;
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
  CLPPSP<WeightT> cli(argc, argv, "point-to-point shortest-path");
  if (!cli.ParseArgs())
    return -1;
  WeightedBuilder b(cli);
  WGraph g = b.MakeGraph();

  dist_array = new WeightT[g.num_nodes()];
  for (int i = 0; i < g.num_nodes(); i++){
    dist_array[i] = kDistInf;
  }

  input_delta = cli.delta();
  pq = new EagerPriorityQueue<WeightT>(dist_array, input_delta);
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
