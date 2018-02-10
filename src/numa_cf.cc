#include <iostream> 
#include <vector>
#include "intrinsics.h"
#include "segmentgraph_cache.hh"

WGraph edges;
typedef double defined_type_0 [ 20]; 
defined_type_0 * __restrict  latent_vec;
typedef double defined_type_1 [ 20]; 
defined_type_1 * __restrict  error_vec;

double step; 
double lambda; 
int K; 
template <typename APPLY_FUNC > VertexSubset<NodeID>* edgeset_apply_pull_parallel_weighted_pull_edge_based_load_balance(WGraph & g , APPLY_FUNC apply_func, GraphSegments<WNode,NodeID>* graphSegments) 
{ 
    long numVertices = g.num_nodes(), numEdges = g.num_edges();
    if (g.offsets_ == nullptr) g.SetUpOffsets(true);
  SGOffset * edge_in_index = g.offsets_;

  /* stage 2 */
  for (int i = 0; i < graphSegments->numSegments; i++) {
    auto sg = graphSegments->getSegmentedGraph(i);
    int* segmentVertexArray = sg->vertexArray;
    WNode* segmentEdgeArray = sg->edgeArray;
#pragma omp parallel for schedule(dynamic, 16)
    for (int localVertexId = 0; localVertexId < sg->numVertices; localVertexId++) {
      int u = sg->graphId[localVertexId];
      int start = segmentVertexArray[localVertexId];
      int end = segmentVertexArray[localVertexId+1];
      for (int neighbor = start; neighbor < end; neighbor++) {
  	WNode v = segmentEdgeArray[neighbor];
  	apply_func ( v.v , u, v.w );
      }
    }
  }

  return new VertexSubset<NodeID>(g.num_nodes(), g.num_nodes());
} //end of edgeset apply function 
void updateEdge(NodeID src, NodeID dst, int rating) 
{
  double estimate = (0) ;
  for ( int i = (0) ; i < K; i++ )
  {
    estimate += (latent_vec[src][i] * latent_vec[dst][i]);
  }
  double err = (rating - estimate);
  for ( int i = (0) ; i < K; i++ )
  {
    error_vec[dst][i] += (latent_vec[src][i] * err);
  }
};
void updateVertex(NodeID v) 
{
  for ( int i = (0) ; i < K; i++ )
  {
    latent_vec[v][i] += (step * (( -lambda * latent_vec[v][i]) + error_vec[v][i]));
    error_vec[v][i] = (0) ;
  }
};
void initVertex(NodeID v) 
{
  for ( int i = (0) ; i < K; i++ )
  {
    latent_vec[v][i] = ((float) 0.5) ;
    error_vec[v][i] = (0) ;
  }
};
int main(int argc, char * argv[] ) 
{
  edges = builtin_loadWeightedEdgesFromFile ( argv[(1) ]) ;
  latent_vec = new defined_type_0 [ builtin_getVertices(edges) ];
  error_vec = new defined_type_1 [ builtin_getVertices(edges) ];
  step = ((float) 3.5e-07) ;
  lambda = ((float) 0.001) ;
  K = (20) ;
  for ( int trail = (0) ; trail < (10) ; trail++ )
  {
    parallel_for (int i = 0; i < builtin_getVertices(edges) ; i++) {
      initVertex(i);
    };
    startTimer() ;

    // Build the segmented graph from graph
    int numSegments = 7;
    int segmentRange = (edges.num_nodes() + numSegments) / numSegments;
    GraphSegments<WNode,int>* graphSegments = new GraphSegments<WNode,int>(numSegments);
    BuildPullSegmentedGraphsWeighted(&edges, graphSegments, segmentRange);

    omp_set_nested(1);
    Timer numa_timer;
    numa_timer.Start();

    for ( int i = (0) ; i < (10) ; i++ )
    {
      edgeset_apply_pull_parallel_weighted_pull_edge_based_load_balance(edges, updateEdge, graphSegments); 
      parallel_for (int i = 0; i < builtin_getVertices(edges) ; i++) {
        updateVertex(i);
      };
#ifdef COUNT
      double latent_sum = 0;
#pragma omp parallel for reduction(+ : latent_sum)
      for (int v = 0; v < builtin_getVertices(edges); v++) {
	for (int i = 0; i < K; i++) {
	  latent_sum += latent_vec[v][i];
	}
      }
      std::cout << "latent_sum=" << latent_sum << endl;
#endif
    }

    numa_timer.Stop();
    PrintTime("NUMA Time", numa_timer.Seconds());

    double elapsed_time = stopTimer() ;
    std::cout << "elapsed time: "<< std::endl;
    std::cout << elapsed_time<< std::endl;
  }
};

