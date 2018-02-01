// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <algorithm>
#include <iostream>
#include <vector>
#include <chrono>
#include <numa.h>

#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "pvector.h"
#include "segmentgraph.hh"

/*
GAP Benchmark Suite
Kernel: PageRank (PR)
Author: Scott Beamer

Will return pagerank scores for all vertices once total change < epsilon

This PR implementation uses the traditional iterative approach. This is done
to ease comparisons to other implementations (often use same algorithm), but
it is not necesarily the fastest way to implement it. It does perform the
updates in the pull direction to remove the need for atomics.
*/


using namespace std;

typedef float ScoreT;
const float kDamp = 0.85;

pvector<ScoreT> PageRankPull(const Graph &g, int max_iters,
                             double epsilon = 0) {
  int num_nodes = g.num_nodes();
  const ScoreT init_score = 1.0f / num_nodes;
  const ScoreT base_score = (1.0f - kDamp) / num_nodes;

  //Build the segmented graph from g
  int num_places = omp_get_num_places();
  int numSegments = num_places * 5;
  int segmentRange = (num_nodes + numSegments) / numSegments;
  GraphSegments<int,int>* graphSegments = new GraphSegments<int,int>(numSegments, num_places);
  BuildPullSegmentedGraphs(&g, graphSegments, segmentRange);

  // Build per socket local buffer
  float *incoming_total[num_places];
  float *outgoing_contrib[num_places];
  for (int socketId = 0; socketId < num_places; socketId++) {
    incoming_total[socketId] = (float *)numa_alloc_onnode(sizeof(float) * num_nodes, socketId);
#pragma omp parallel for
    for (int n = 0; n < num_nodes; n++) {
      incoming_total[socketId][n] = 0;
    }

    outgoing_contrib[socketId] = (float *)numa_alloc_onnode(sizeof(float) * num_nodes, socketId);
  }

#ifdef LOAD_MSG
  for (int segmentId = 0; segmentId < numSegments; segmentId++) {
    auto sg = graphSegments->getSegmentedGraph(segmentId);
    cout << "segmentId=" << segmentId <<  " numVertices=" << graphSegments->getSegmentedGraph(segmentId)->numVertices
	 << " numEdges=" << graphSegments->getSegmentedGraph(segmentId)->numEdges << endl;
  }
#endif

  Timer numa_timer;
  numa_timer.Start();

  pvector<ScoreT> scores(num_nodes, init_score);

  for (int iter=0; iter < max_iters; iter++) {
    double error = 0;

#ifdef TIME_MSG
    Timer debug_timer;
    debug_timer.Start();
#endif
    /* stage 1: compute outgoing_contrib, global sequential */
//     for (int segmentId = 0; segmentId < numSegments; segmentId++) {
//       int socketId = segmentId % num_places;
// #pragma omp parallel for
//       for (int i = 0; i < segmentRange; i++) {
// 	int n = segmentId * segmentRange + i;
// 	if (n < num_nodes)
// 	  outgoing_contrib[socketId][n] = scores[n] / g.out_degree(n);
//       }
//     }
#pragma omp parallel for
    for (int n = 0; n < num_nodes; n++) {
      for (int i = 0; i < num_places; i++) {
	outgoing_contrib[i][n] = scores[n] / g.out_degree(n);
      }
    }
    graphSegments->distribute();
#ifdef TIME_MSG
    debug_timer.Stop();
    cout << "stage 1 took " << debug_timer.Seconds() << " seconds" << endl;
    debug_timer.Start();
#endif

    /* stage 2: pull from neighbor, used to be global random, now local random */
    omp_set_nested(1);
#ifdef TIME_MSG
  auto start = chrono::steady_clock::now();
#endif
#pragma omp parallel num_threads(num_places) proc_bind(spread)
    {
      int socketId = omp_get_place_num();
      int n_procs = omp_get_place_num_procs(socketId);
      while (true) {
	auto sg = graphSegments->getSegmentedGraphFromGroup(socketId);
	if (!sg)
	  break;
	int* segmentVertexArray = sg->vertexArray;
	int* segmentEdgeArray = sg->edgeArray;
#pragma omp parallel num_threads(n_procs) proc_bind(close)
	{
#pragma omp for schedule(dynamic, 64)
	  for (int localVertexId = 0; localVertexId < sg->numVertices; localVertexId++) {
	    int u = sg->graphId[localVertexId];
	    int start = segmentVertexArray[localVertexId];
	    int end = segmentVertexArray[localVertexId+1];
	    for (int neighbor = start; neighbor < end; neighbor++) {
	      int v = segmentEdgeArray[neighbor];
	      incoming_total[socketId][u] += outgoing_contrib[socketId][v];
	    }
	  }
	}
      }
#ifdef TIME_MSG
#pragma omp critical
      {
	auto end = chrono::steady_clock::now();
	cout << "socket=" << socketId << " time=" << chrono::duration_cast<chrono::duration<double>>(end - start).count()
	     << "seconds" << endl;
      }
#endif
    }
#ifdef TIME_MSG
    debug_timer.Stop();
    cout << "stage 2 took " << debug_timer.Seconds() << " seconds" << endl;
    debug_timer.Start();
#endif

    /* stage 3 reduce scores, global sequential */
#pragma omp parallel for reduction(+ : error)
    for (NodeID n=0; n < num_nodes; n++) {
      ScoreT old_score = scores[n];
      ScoreT global_incoming_total = 0;
      for (int socketId = 0; socketId < num_places; socketId++) {
	global_incoming_total += incoming_total[socketId][n];
	incoming_total[socketId][n] = 0;
      }
      scores[n] = base_score + kDamp * global_incoming_total;
      error += fabs(scores[n] - old_score);
    }

#ifdef TIME_MSG
    debug_timer.Stop();
    cout << "stage 3 took " << debug_timer.Seconds() << " seconds" << endl;
#endif

    printf(" %2d    %lf\n", iter, error);
    if (error < epsilon)
      break;
  }

  numa_timer.Stop();
  PrintTime("NUMA Time", numa_timer.Seconds());

  for (int socketId = 0; socketId < num_places; socketId++) {
    numa_free(incoming_total[socketId], sizeof(float) * num_nodes);
    numa_free(outgoing_contrib[socketId], sizeof(float) * num_nodes);
  }

  return scores;
}


void PrintTopScores(const Graph &g, const pvector<ScoreT> &scores) {
  vector<pair<NodeID, ScoreT>> score_pairs(g.num_nodes());
  for (NodeID n=0; n < g.num_nodes(); n++) {
    score_pairs[n] = make_pair(n, scores[n]);
  }
  int k = 5;
  vector<pair<ScoreT, NodeID>> top_k = TopK(score_pairs, k);
  k = min(k, static_cast<int>(top_k.size()));
  for (auto kvp : top_k)
    cout << kvp.second << ":" << kvp.first << endl;
}


// Verifies by asserting a single serial iteration in push direction has
//   error < target_error
bool PRVerifier(const Graph &g, const pvector<ScoreT> &scores,
                        double target_error) {
  const ScoreT base_score = (1.0f - kDamp) / g.num_nodes();
  pvector<ScoreT> incomming_sums(g.num_nodes(), 0);
  double error = 0;

    float total_scores = 0;
    for (NodeID u : g.vertices()){
        total_scores += scores[u];
    }
    std::cout << "total scores: " << total_scores << std::endl;

  for (NodeID u : g.vertices()) {
    ScoreT outgoing_contrib = scores[u] / g.out_degree(u);
    for (NodeID v : g.out_neigh(u))
      incomming_sums[v] += outgoing_contrib;
  }
  for (NodeID n : g.vertices()) {
    error += fabs(base_score + kDamp * incomming_sums[n] - scores[n]);
    incomming_sums[n] = 0;
  }
  PrintTime("Total Error", error);
  return error < target_error;
}


int main(int argc, char* argv[]) {
  CLPageRank cli(argc, argv, "pagerank", 1e-4, 20);
  if (!cli.ParseArgs())
    return -1;
  Builder b(cli);
  Graph g = b.MakeGraph();
  auto PRBound = [&cli] (const Graph &g) {
    return PageRankPull(g, cli.max_iters(), cli.tolerance());
  };

  auto VerifierBound = [&cli] (const Graph &g, const pvector<ScoreT> &scores) {
    return PRVerifier(g, scores, cli.tolerance());
  };
  BenchmarkKernel(cli, g, PRBound, PrintTopScores, VerifierBound);
  return 0;
}
