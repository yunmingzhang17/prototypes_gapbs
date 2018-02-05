#include <iostream> 
#include <vector>

#include "intrinsics.h"
#include "segmentgraph.hh"

Graph edges; 
typedef struct struct_delta_out_degree { 
  double delta;
  int out_degree;
} struct_delta_out_degree;

static double **local_ngh_sum;
static struct_delta_out_degree **local_array_of_struct_delta_out_degree;

double  * __restrict cur_rank;
double  * __restrict ngh_sum;
struct_delta_out_degree  * __restrict array_of_struct_delta_out_degree;
int  * __restrict generated_tmp_vector_3;
double damp; 
double beta_score; 
double epsilon2; 
double epsilon; 
template <typename APPLY_FUNC , typename PUSH_APPLY_FUNC> VertexSubset<NodeID>* edgeset_apply_hybrid_dense_parallel_from_vertexset_pull_frontier_bitvector(Graph & g , VertexSubset<NodeID>* from_vertexset, APPLY_FUNC apply_func, PUSH_APPLY_FUNC push_apply_func, GraphSegments<int, int>* graphSegments) 
{ 
    long numVertices = g.num_nodes(), numEdges = g.num_edges();
    from_vertexset->toSparse();
    long m = from_vertexset->size();
    // used to generate nonzero indices to get degrees
    uintT *degrees = newA(uintT, m);
    // We probably need this when we get something that doesn't have a dense set, not sure
    // We can also write our own, the eixsting one doesn't quite work for bitvectors
    //from_vertexset->toSparse();
    {
        parallel_for (long i = 0; i < m; i++) {
            NodeID v = from_vertexset->dense_vertex_set_[i];
            degrees[i] = g.out_degree(v);
        }
    }
    uintT outDegrees = sequence::plusReduce(degrees, m);
    if (m + outDegrees > numEdges / 20) {
    from_vertexset->toDense();
  Bitmap bitmap(numVertices);
  bitmap.reset();
  parallel_for(int i = 0; i < numVertices; i+=32){
     int start = i;
     int end = (((i + 32) < numVertices)? (i+32):numVertices);
     for(int j = start; j < end; j++){
        if (from_vertexset->bool_map_[j])
          bitmap.set_bit(j);
     }
  }

    /* stage 1: compute outgoing_contrib, global sequential */
#ifdef TIME_MSG
    Timer debug_timer;
    debug_timer.Start();
#endif
   int num_places = omp_get_num_places();
#pragma omp parallel for
   for (int n = 0; n < g.num_nodes(); n++) {
      for (int i = 0; i < num_places; i++) {
        local_array_of_struct_delta_out_degree[i][n].delta = array_of_struct_delta_out_degree[n].delta;
        local_array_of_struct_delta_out_degree[i][n].out_degree = array_of_struct_delta_out_degree[n].out_degree;
      }
    }
    graphSegments->resetGroups();
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
	      if (bitmap.get_bit(v)) { 
		apply_func(v, u, socketId);
	      }
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

// #pragma omp parallel for schedule(dynamic, 64)
//     for ( NodeID d=0; d < g.num_nodes(); d++) {
//       for(NodeID s : g.in_neigh(d)){
//         if (bitmap.get_bit(s)) { 
//           apply_func ( s , d  );
//         }
//       } //end of loop on in neighbors
//     } //end of outer for loop

    /* stage 3 reduce scores, global sequential */
#pragma omp parallel for
    for (NodeID n=0; n < g.num_nodes(); n++) {
      double global_ngh_sum = 0;
      for (int socketId = 0; socketId < num_places; socketId++) {
	global_ngh_sum += local_ngh_sum[socketId][n];
	local_ngh_sum[socketId][n] = 0;
      }
      ngh_sum[n] = global_ngh_sum;
    }

#ifdef TIME_MSG
    debug_timer.Stop();
    cout << "stage 3 took " << debug_timer.Seconds() << " seconds" << endl;
#endif

    return new VertexSubset<NodeID>(g.num_nodes(), g.num_nodes());
} else {
      parallel_for (long i=0; i < m; i++) {
    NodeID s = from_vertexset->dense_vertex_set_[i];
    int j = 0;
        for(NodeID d : g.out_neigh(s)){
          push_apply_func ( s , d  );
        } //end of for loop on neighbors
      }
      return new VertexSubset<NodeID>(g.num_nodes(), g.num_nodes());
} //end of else
} //end of edgeset apply function 
void updateEdge_push_ver(NodeID src, NodeID dst) 
{
  writeAdd( &ngh_sum[dst], (array_of_struct_delta_out_degree[src].delta  / array_of_struct_delta_out_degree[src].out_degree ) ); 
};
void generated_vector_op_apply_func_4(NodeID v) 
{
  array_of_struct_delta_out_degree[v].out_degree  = generated_tmp_vector_3[v];
};
void delta_generated_vector_op_apply_func_2(NodeID v) 
{
  array_of_struct_delta_out_degree[v].delta  = (((float) 1)  / builtin_getVertices(edges) );
};
void ngh_sum_generated_vector_op_apply_func_1(NodeID v) 
{
  ngh_sum[v] = ((float) 0) ;
};
void cur_rank_generated_vector_op_apply_func_0(NodeID v) 
{
  cur_rank[v] = (0) ;
};
void updateEdge(NodeID src, NodeID dst, int socketId)
{
  local_ngh_sum[socketId][dst] += local_array_of_struct_delta_out_degree[socketId][src].delta / local_array_of_struct_delta_out_degree[socketId][src].out_degree;
  //ngh_sum[dst] += (array_of_struct_delta_out_degree[src].delta  / array_of_struct_delta_out_degree[src].out_degree );
};
bool updateVertexFirstRound(NodeID v) 
{
  bool output ;
  array_of_struct_delta_out_degree[v].delta  = ((damp * ngh_sum[v]) + beta_score);
  cur_rank[v] += array_of_struct_delta_out_degree[v].delta ;
  array_of_struct_delta_out_degree[v].delta  = (array_of_struct_delta_out_degree[v].delta  - (((float) 1)  / builtin_getVertices(edges) ));
  output = (fabs(array_of_struct_delta_out_degree[v].delta ) ) > ((epsilon2 * cur_rank[v]));
  ngh_sum[v] = (0) ;
  return output;
};
bool updateVertex(NodeID v) 
{
  bool output ;
  array_of_struct_delta_out_degree[v].delta  = (ngh_sum[v] * damp);
  cur_rank[v] += array_of_struct_delta_out_degree[v].delta ;
  ngh_sum[v] = (0) ;
  output = (fabs(array_of_struct_delta_out_degree[v].delta ) ) > ((epsilon2 * cur_rank[v]));
  return output;
};
void printRank(NodeID v) 
{
  std::cout << cur_rank[v]<< std::endl;
};
int main(int argc, char * argv[] ) 
{
  edges = builtin_loadEdgesFromFile ( argv[(1) ]) ;
  cur_rank = new double [ builtin_getVertices(edges) ];
  ngh_sum = new double [ builtin_getVertices(edges) ];
  array_of_struct_delta_out_degree = new struct_delta_out_degree [ builtin_getVertices(edges) ];
  generated_tmp_vector_3 = new int [ builtin_getVertices(edges) ];
  damp = ((float) 0.85) ;
  beta_score = ((((float) 1)  - damp) / builtin_getVertices(edges) );
  epsilon2 = ((float) 0.1) ;
  epsilon = ((float) 1e-07) ;
  parallel_for (int i = 0; i < builtin_getVertices(edges) ; i++) {
    cur_rank_generated_vector_op_apply_func_0(i);
  };
  parallel_for (int i = 0; i < builtin_getVertices(edges) ; i++) {
    ngh_sum_generated_vector_op_apply_func_1(i);
  };
  parallel_for (int i = 0; i < builtin_getVertices(edges) ; i++) {
    delta_generated_vector_op_apply_func_2(i);
  };
  generated_tmp_vector_3 = builtin_getOutDegrees(edges) ;
  parallel_for (int i = 0; i < builtin_getVertices(edges) ; i++) {
    generated_vector_op_apply_func_4(i);
  };
  int n = builtin_getVertices(edges) ;
  VertexSubset<int> *  frontier = new VertexSubset<int> ( builtin_getVertices(edges)  , n);
  startTimer() ;

  // Build the segmented graph from edges
  int num_nodes = edges.num_nodes();
  int num_places = omp_get_num_places();
  int numSegments = num_places * 5;
  int segmentRange = (num_nodes + numSegments) / numSegments;
  GraphSegments<int,int>* graphSegments = new GraphSegments<int,int>(numSegments, num_places);
  BuildPullSegmentedGraphs(&edges, graphSegments, segmentRange);

  // Build per socket local buffer
  local_ngh_sum = new double*[num_places];
  local_array_of_struct_delta_out_degree = new struct_delta_out_degree*[num_places];
  for (int socketId = 0; socketId < num_places; socketId++) {
    local_ngh_sum[socketId] = (double *)numa_alloc_onnode(sizeof(double) * num_nodes, socketId);
#pragma omp parallel for
    for (int n = 0; n < num_nodes; n++) {
      local_ngh_sum[socketId][n] = 0;
    }
    local_array_of_struct_delta_out_degree[socketId] = (struct_delta_out_degree *)numa_alloc_onnode(sizeof(struct_delta_out_degree) * num_nodes, socketId);
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

  for ( int i = (1) ; i < (11) ; i++ )
  {
  edgeset_apply_hybrid_dense_parallel_from_vertexset_pull_frontier_bitvector(edges, frontier, updateEdge, updateEdge_push_ver,
    graphSegments); 
    if ((i) == ((1) ))
     { 
      auto ____graphit_tmp_out = new VertexSubset <NodeID> ( builtin_getVertices(edges)  , 0 );
bool * next5 = newA(bool, builtin_getVertices(edges) );
      parallel_for (int v = 0; v < builtin_getVertices(edges) ; v++) {
        next5[v] = 0;
if ( updateVertexFirstRound( v ) )
          next5[v] = 1;
      } //end of loop
____graphit_tmp_out->num_vertices_ = sequence::sum( next5, builtin_getVertices(edges)  );
____graphit_tmp_out->bool_map_ = next5;

      frontier  = ____graphit_tmp_out; 
     } 
    else
     { 
      auto ____graphit_tmp_out = new VertexSubset <NodeID> ( builtin_getVertices(edges)  , 0 );
bool * next6 = newA(bool, builtin_getVertices(edges) );
      parallel_for (int v = 0; v < builtin_getVertices(edges) ; v++) {
        next6[v] = 0;
if ( updateVertex( v ) )
          next6[v] = 1;
      } //end of loop
____graphit_tmp_out->num_vertices_ = sequence::sum( next6, builtin_getVertices(edges)  );
____graphit_tmp_out->bool_map_ = next6;

      frontier  = ____graphit_tmp_out; 

     } 
    std::cout << "frontier size=" << builtin_getVertexSetSize(frontier) << std::endl;
  }

  numa_timer.Stop();
  PrintTime("NUMA Time", numa_timer.Seconds());

  for (int socketId = 0; socketId < num_places; socketId++) {
    numa_free(local_ngh_sum[socketId], sizeof(double) * num_nodes);
    numa_free(local_array_of_struct_delta_out_degree[socketId], sizeof(struct_delta_out_degree) * num_nodes);
  }

  double elapsed_time = stopTimer() ;
  std::cout << "elapsed time: "<< std::endl;
  std::cout << elapsed_time<< std::endl;
};

