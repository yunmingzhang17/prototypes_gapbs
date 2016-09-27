// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <vector>

#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "pvector.h"
#include "segmentgraph.hh"

/*
GAP Benchmark Suite
Kernel: Triangle Counting (TC)
Author: Scott Beamer

Will count the number of triangles (cliques of size 3)

Requires input graph:
  - to be undirected
  - no duplicate edges (or else will be counted as multiple triangles)
  - neighborhoods are sorted by vertex identifiers

Other than symmetrizing, the rest of the requirements are done by SquishCSR
during graph building.

This implementation reduces the search space by counting each triangle only
once. A naive implementation will count the same triangle six times because
each of the three vertices (u, v, w) will count it in both ways. To count
a triangle only once, this implementation only counts a triangle if u > v > w.
Once the remaining unexamined neighbors identifiers get too big, it can break
out of the loop, but this requires that the neighbors to be sorted.

Another optimization this implementation has is to relabel the vertices by
degree. This is beneficial if the average degree is high enough and if the
degree distribution is sufficiently non-uniform. To decide whether or not
to relabel the graph, we use the heuristic in WorthRelabelling.
*/


using namespace std;

size_t OrderedCount(const Graph &g) {
  
  size_t total = 0;
  /*
  //#pragma omp parallel for reduction(+ : total) schedule(dynamic, 64)
  for (NodeID u=0; u < g.num_nodes(); u++) {
    for (NodeID v : g.out_neigh(u)) {
#ifdef DEBUG
      cout << "edge from u: " << u << " to v: " << v << endl;
#endif
      if (v > u)
        break;
      auto it = g.out_neigh(u).begin();
      for (NodeID w : g.out_neigh(v)) {
        if (w > v)
          break;
        while (*it < w)
          it++;
        if (w == *it){
#ifdef DEBUG
      cout << "triangle: u, v, w: " << u << " " << v << " "  << w << endl;
#endif
	  total++;
        }
      }
    }
  }
  
  total = 0;
  */

  //Build the segmented graph from g

  //4 M integers for 32 M LLC cache, using about 50% of the cache
  int numElementsPerSegment = 1024*1024*4; 
  //4 k integers for 32 K L1 cache, usign about 50% of the cache
  int numElementsPerSlice = 1024*4;
  int numSegments = (g.num_edges() + numElementsPerSegment)/numElementsPerSegment; 
  int numSlices = g.num_nodes()/numElementsPerSlice;

  cout << "number of segments: " << numSegments << endl;

  GraphSegments<int,int>* graphSegments = new GraphSegments<int,int>(numSegments, numSlices, numElementsPerSlice);

  BuildCacheSegmentedGraphs(&g, graphSegments, numSlices, numElementsPerSegment, numElementsPerSlice);

  //Perform computation within each segment
  for (int segmentId = 0; segmentId < numSegments; segmentId++){
    auto sg = graphSegments->getSegmentedGraph(segmentId);
    int* segmentVertexArray = sg->vertexArray;
    int* segmentEdgeArray = sg->edgeArray;
    size_t local_total = 0;
    int localVertexId;

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
#endif

    #pragma omp parallel for reduction(+ : local_total) schedule(dynamic, 64)   
    for (localVertexId = 0; localVertexId < sg->numVertices; localVertexId++){
      int u = sg->graphId[localVertexId];
      int start = segmentVertexArray[localVertexId];
      int end = segmentVertexArray[localVertexId+1];
      for (int neighbor = start; neighbor < end; neighbor++){
	int v = segmentEdgeArray[neighbor];
#ifdef DEBUG2
	cout << "edge from u: " << u << " to v: " << v << endl;
#endif
	if (v > u)
	  break;
	auto it = g.out_neigh(u).begin();                                     
	for (NodeID w : g.out_neigh(v)) {                                     
	  if (w > v)                                                          
	    break;                                                            
	  while (*it < w)                                                     
	    it++;                                                             
	  if (w == *it){                                                      
#ifdef DEBUG2
      cout << "triangle: u, v, w: " << u << " " <<  v << " "  << w << endl;
#endif
	    local_total++;                                                    
           }
        }                                                                     
      }
      
    }//end of localVertexId for
    total += local_total;
  }// end of segment for

  return total;

}


// heuristic to see if sufficently dense power-law graph
bool WorthRelabelling(const Graph &g) {
  int64_t average_degree = g.num_edges() / g.num_nodes();
  if (average_degree < 10)
    return false;
  SourcePicker<Graph> sp(g);
  int64_t num_samples = min(int64_t(1000), g.num_nodes());
  int64_t sample_total = 0;
  pvector<int64_t> samples(num_samples);
  for (int64_t trial=0; trial < num_samples; trial++) {
    samples[trial] = g.out_degree(sp.PickNext());
    sample_total += samples[trial];
  }
  sort(samples.begin(), samples.end());
  double sample_average = static_cast<double>(sample_total) / num_samples;
  double sample_median = samples[num_samples/2];
  return sample_average / 2 > sample_median;
}


// uses heuristic to see if worth relabeling
size_t Hybrid(const Graph &g) {
  if (WorthRelabelling(g))
    return OrderedCount(Builder::RelabelByDegree(g));
  else
    return OrderedCount(g);
}


void PrintTriangleStats(const Graph &g, size_t total_triangles) {
  cout << total_triangles << " triangles" << endl;
}


// Compares with simple serial implementation that uses std::set_intersection
bool TCVerifier(const Graph &g, size_t test_total) {
  size_t total = 0;
  vector<NodeID> intersection;
  intersection.reserve(g.num_nodes());
  for (NodeID u : g.vertices()) {
    for (NodeID v : g.out_neigh(u)) {
      auto new_end = set_intersection(g.out_neigh(u).begin(),
                                      g.out_neigh(u).end(),
                                      g.out_neigh(v).begin(),
                                      g.out_neigh(v).end(),
                                      intersection.begin());
      intersection.resize(new_end - intersection.begin());
      total += intersection.size();
    }
  }
  total = total / 6;  // each triangle was counted 6 times
  if (total != test_total)
    cout << total << " != " << test_total << endl;
  return total == test_total;
}


int main(int argc, char* argv[]) {
  CLApp cli(argc, argv, "triangle count");
  if (!cli.ParseArgs())
    return -1;
  Builder b(cli);
  Graph g = b.MakeGraph();
  if (g.directed()) {
    cout << "Input graph is directed but tc requires undirected" << endl;
    return -2;
  }
  BenchmarkKernel(cli, g, Hybrid, PrintTriangleStats, TCVerifier);
  return 0;
}
