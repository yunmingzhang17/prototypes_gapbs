// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <algorithm>
#include <iostream>
#include <vector>
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
  const ScoreT init_score = 1.0f / g.num_nodes();
  const ScoreT base_score = (1.0f - kDamp) / g.num_nodes();

  //Build the segmented graph from g
  //4 M integers for 32 M LLC cache, using about 50% of the cache
  int numElementsPerSegment = 1024*1024*4; 
  //4 k integers for 32 K L1 cache, usign about 50% of the cache
  int numElementsPerSlice = 1024*4;
  int numSegments = (2*g.num_edges() + numElementsPerSegment)/numElementsPerSegment; 
  int numSlices = g.num_nodes()/numElementsPerSlice;
  cout << "number of segments: " << numSegments << endl;
  GraphSegments<int,int>* graphSegments = new GraphSegments<int,int>(numSegments, numSlices, numElementsPerSlice);
  BuildCacheSegmentedGraphs(&g, graphSegments, numSlices, numElementsPerSegment, numElementsPerSlice);

  pvector<ScoreT> scores(g.num_nodes(), init_score);
  //pvector<ScoreT> outgoing_contrib(g.num_nodes());

  for (int iter=0; iter < max_iters; iter++) {
    double error = 0;

    //omp_set_nested(1);
    //#pragma omp parallel for
    //Perform computation within each segment
    for (int segmentId = 0; segmentId < numSegments; segmentId++) {
      auto sg = graphSegments->getSegmentedGraph(segmentId);
      int* segmentVertexArray = sg->vertexArray;
      int* segmentEdgeArray = sg->edgeArray;
      int localVertexId;

      nodemask_t mask;
      struct bitmask *bmaskp = numa_bitmask_alloc(num_numa_node);
      nodemask_zero(&mask);
      nodemask_set_compat(&mask, segmentId % num_numa_node);
      copy_nodemask_to_bitmask(&mask, bmaskp);
      numa_bind(bmaskp);

#ifdef DEBUG2
      cout << "in segment: " << segmentId << endl;
      cout << "num vertices: " << sg->numVertices << " num edges: " << sg->numEdges << endl;
      cout << "segment vertex array: " << endl;
      for (int i = 0; i < sg->numVertices + 1; i++){
	cout << " " << segmentVertexArray[i];
      }
      cout << endl;
      cout << "segment edge array: " << endl;
      for (int i = 0; i < sg->numEdges; i++){
	cout << " " << segmentEdgeArray[i];
      }
      cout << endl;
      cout << "running on node mask: " << mask << endl;
#endif

#pragma omp parallel for
      for (localVertexId = 0; localVertexId < sg->numVertices; localVertexId++) {
	int n = sg->graphId[localVertexId];
	sg->outgoing_contrib[n] = sg->scores[n] / g.out_degree(n);
      }
      
#pragma omp parallel for reduction(+ : error) schedule(dynamic, 64)
      for (localVertexId = 0; localVertexId < sg->numVertices; localVertexId++) {
	int u = sg->graphId[localVertexId];
	ScoreT incoming_total = 0;
	int start = segmentVertexArray[localVertexId];
	int end = segmentVertexArray[localVertexId+1];
	for (int neighbor = start; neighbor < end; neighbor++) {
	  int v = segmentEdgeArray[neighbor];
#ifdef DEBUG2
	  cout << "edge from v: " << v << " to u: " << u << endl;
#endif
	  incoming_total += sg->outgoing_contrib[v];
	}
	ScoreT old_score = sg->scores[u];
	sg->scores[u] = base_score + kDamp * incoming_total;
	error += fabs(sg->scores[u] - old_score);
      }
    }
    printf(" %2d    %lf\n", iter, error);
    if (error < epsilon)
      break;
  }
  
  // copy scores to the final array
  for (int segmentId = 0; segmentId < numSegments; segmentId++) {
    auto sg = graphSegments->getSegmentedGraph(segmentId);
#pragma omp parallel for
    for (int localVertexId = 0; localVertexId < sg->numVertices; localVertexId++) {
      int n = sg->graphId[localVertexId];
      scores[n] = sg->scores[n];
    }
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
