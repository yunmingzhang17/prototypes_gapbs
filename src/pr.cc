// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <algorithm>
#include <iostream>
#include <vector>
#include <numa.h>
#include <omp.h>
#include <assert.h>

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

static int threads;
static int sockets;
static int threads_per_socket;

void numa_thread_init(void) {
  threads = numa_num_configured_cpus();
  sockets = numa_num_configured_nodes();
  threads_per_socket = threads / sockets;
  assert(numa_available() != -1);
  omp_set_dynamic(0);
  omp_set_num_threads(threads);

#pragma omp parallel
  {
    int thread_id = omp_get_thread_num();
    int socket_id = thread_id / threads_per_socket;
    assert(numa_run_on_node(socket_id) == 0);

#ifdef DEBUG_MSG
    cout << "binding: omp_get_num_threads=" << omp_get_num_threads() << endl;
#endif
  }
}

pvector<ScoreT> PageRankPull(const Graph &g, int max_iters,
                             double epsilon = 0) {
  const ScoreT init_score = 1.0f / g.num_nodes();
  const ScoreT base_score = (1.0f - kDamp) / g.num_nodes();
  pvector<ScoreT> scores(g.num_nodes(), init_score);
  pvector<ScoreT> outgoing_contrib(g.num_nodes());

  /* bind memory */
  int num_nodes = g.num_nodes();
  int numSegments = num_numa_node;
  int segmentRange = (num_nodes + numSegments) / numSegments;
  GraphSegments<int,int>* graphSegments = new GraphSegments<int,int>(numSegments, segmentRange, num_nodes);
  BuildCacheSegmentedGraphs(&g, graphSegments, segmentRange);

#ifdef LOAD_MSG
  cout << "socket 0 has " << graphSegments->getSegmentedGraph(0)->numDstVertices << " vertices" << endl;
  cout << "socket 0 has " << graphSegments->getSegmentedGraph(0)->numEdges << " edges" << endl;
  cout << "socket 1 has " << graphSegments->getSegmentedGraph(1)->numDstVertices << " vertices" << endl;
  cout << "socket 1 has " << graphSegments->getSegmentedGraph(1)->numEdges << " edges" << endl;
#endif

  /* bind threads */
  numa_thread_init();

  Timer numa_timer;
  numa_timer.Start();

  for (int iter=0; iter < max_iters; iter++) {
    double error = 0;


#ifdef TIME_MSG
    Timer debug_timer;
    debug_timer.Start();
#endif

    /* stage 1: compute outgoing_contrib, global sequential */
#pragma omp parallel for
    for (NodeID n=0; n < g.num_nodes(); n++)
      outgoing_contrib[n] = scores[n] / g.out_degree(n);

#ifdef TIME_MSG
    debug_timer.Stop();
    cout << "stage 1 took " << debug_timer.Seconds() << " seconds" << endl;
    debug_timer.Start();
#endif

    /* stage 2: pull from neighbor, used to be global random, now local random */
#pragma omp parallel
  {
    int thread_id = omp_get_thread_num();
    int socket_id = thread_id / threads_per_socket;
    auto sg = graphSegments->getSegmentedGraph(socket_id);
    int* segmentVertexArray = sg->dstVertexArray;
    int* segmentEdgeArray = sg->edgeArray;

    int vertices_per_thread = (sg->numDstVertices + threads_per_socket) / threads_per_socket;
    int start_index = (thread_id % threads_per_socket) * vertices_per_thread;
    int end_index = min(start_index + vertices_per_thread, sg->numDstVertices);

    for (int localVertexId = start_index; localVertexId < end_index; localVertexId++) {
      int u = sg->graphId[localVertexId];
      int start_neighbor = segmentVertexArray[localVertexId];
      int end_neighbor = segmentVertexArray[localVertexId+1];
      for (int neighbor = start_neighbor; neighbor < end_neighbor; neighbor++) {
	int v = segmentEdgeArray[neighbor];
	sg->incoming_total[u] += outgoing_contrib[v];
      }
    }
#ifdef DEBUG_MSG
    cout << "stage 2: omp_get_num_threads=" << omp_get_num_threads() << endl;
#endif
  }

#ifdef TIME_MSG
    debug_timer.Stop();
    cout << "stage 2 took " << debug_timer.Seconds() << " seconds" << endl;
    debug_timer.Start();
#endif

  /* stage 3 reduce scores, global sequential */
#pragma omp parallel for reduction(+ : error) schedule(dynamic, 64)
  for (NodeID n=0; n < num_nodes; n++) {
    ScoreT old_score = scores[n];
    ScoreT global_incoming_total = 0;
    for (int segmentId = 0; segmentId < numSegments; segmentId++) {
      auto sg = graphSegments->getSegmentedGraph(segmentId);
      global_incoming_total += sg->incoming_total[n];
      sg->incoming_total[n] = 0;
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
