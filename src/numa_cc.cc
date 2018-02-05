#include <iostream> 
#include <vector>
#include <unordered_map>
#include <chrono>

#include "intrinsics.h"
#include "segmentgraph.hh"

Graph edges; 
int  * __restrict IDs;
static int **localIDs;

template <typename APPLY_FUNC , typename PUSH_APPLY_FUNC> VertexSubset<NodeID>* edgeset_apply_hybrid_dense_parallel_deduplicatied_from_vertexset_with_frontier(Graph & g , VertexSubset<NodeID>* from_vertexset, APPLY_FUNC apply_func, PUSH_APPLY_FUNC push_apply_func, GraphSegments<int,int>* graphSegments)
{ 
    long numVertices = g.num_nodes(), numEdges = g.num_edges();
    long m = from_vertexset->size();
    from_vertexset->toSparse();
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
  VertexSubset<NodeID> *next_frontier = new VertexSubset<NodeID>(g.num_nodes(), 0);
  bool * next = newA(bool, g.num_nodes());
  parallel_for (int i = 0; i < numVertices; i++)next[i] = 0;
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

  /* stage 1: broadcast global IDs to per socket localIDs */
#ifdef TIME_MSG
  Timer debug_timer;
  debug_timer.Start();
#endif
   int num_places = omp_get_num_places();
#pragma omp parallel for
  for (NodeID n = 0; n < numVertices; n++) {
    for (int i = 0; i < num_places; i++) {
      localIDs[i][n] = IDs[n];
    }
  }
  graphSegments->resetGroups();
#ifdef TIME_MSG
  debug_timer.Stop();
  cout << "stage 1 took " << debug_timer.Seconds() << " seconds" << endl;
  debug_timer.Start();
#endif

  /* stage 2: pull local reads and writes of sg->IDs */
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
	      if (apply_func(v, u, socketId)) {
		next[u] = 1;
	      }
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
    
  /* stage 3: merge localIDs to the global IDs */
#pragma omp parallel for
  for (NodeID n=0; n < g.num_nodes(); n++) {
    int newID = n;
    for (int socketId = 0; socketId < num_places; socketId++) {
      newID = min(newID, localIDs[socketId][n]);
    }
    IDs[n] = newID;
  }
#ifdef TIME_MSG
  debug_timer.Stop();
  cout << "stage 3 took " << debug_timer.Seconds() << " seconds" << endl << endl;;
#endif

  next_frontier->num_vertices_ = sequence::sum(next, numVertices);
  next_frontier->bool_map_ = next;
  return next_frontier;
} else {
    if (g.flags_ == nullptr){
      g.flags_ = new int[numVertices]();
    }
  parallel_for(int i = 0; i < m; i++){
     g.flags_[from_vertexset->dense_vertex_set_[i]] = 0;
  }
    VertexSubset<NodeID> *next_frontier = new VertexSubset<NodeID>(g.num_nodes(), 0);
    if (numVertices != from_vertexset->getVerticesRange()) {
        cout << "edgeMap: Sizes Don't match" << endl;
        abort();
    }
    if (outDegrees == 0) return next_frontier;
    uintT *offsets = degrees;
    long outEdgeCount = sequence::plusScan(offsets, degrees, m);
    uintE *outEdges = newA(uintE, outEdgeCount);
      parallel_for (long i=0; i < m; i++) {
    NodeID s = from_vertexset->dense_vertex_set_[i];
    uintT offset = offsets[i];
    int j = 0;
        for(NodeID d : g.out_neigh(s)){
          if( push_apply_func ( s , d  ) && CAS(&(g.flags_[d]), 0, 1)  ) { 
            outEdges[offset + j] = d; 
          } else { outEdges[offset + j] = UINT_E_MAX; }
          j++;
        } //end of for loop on neighbors
      }
  uintE *nextIndices = newA(uintE, outEdgeCount);
  long nextM = sequence::filter(outEdges, nextIndices, outEdgeCount, nonMaxF());
  free(outEdges);
  free(degrees);
  next_frontier->num_vertices_ = nextM;
  next_frontier->dense_vertex_set_ = nextIndices;
  return next_frontier;
} //end of else
} //end of edgeset apply function 
bool updateEdge_push_ver(NodeID src, NodeID dst) 
{
  bool output4 ;
  bool IDs_trackving_var_3 = (bool) 0;
  IDs_trackving_var_3 = writeMin( &IDs[dst], IDs[src] ); 
  output4 = IDs_trackving_var_3;
  return output4;
};
void IDs_generated_vector_op_apply_func_0(NodeID v) 
{
  IDs[v] = (1) ;
};
bool updateEdge(NodeID src, NodeID dst, int socketId)
{
  bool output2 ;
  bool IDs_trackving_var_1 = (bool) 0;
  if (localIDs[socketId][dst] > localIDs[socketId][src]) {
    localIDs[socketId][dst] = localIDs[socketId][src];
    IDs_trackving_var_1 = true;
  }
  output2 = IDs_trackving_var_1;
  return output2;
};
void init(NodeID v) 
{
  IDs[v] = v;
};
int main(int argc, char * argv[] ) 
{
  edges = builtin_loadEdgesFromFile ( argv[(1) ]) ;
  IDs = new int [ builtin_getVertices(edges) ];
  parallel_for (int i = 0; i < builtin_getVertices(edges) ; i++) {
    IDs_generated_vector_op_apply_func_0(i);
  };
  
  int n = builtin_getVertices(edges) ;
  for ( int trail = (0) ; trail < (10) ; trail++ )
  {
    startTimer() ;
    VertexSubset<int> *  frontier = new VertexSubset<int> ( builtin_getVertices(edges)  , n);
    parallel_for (int i = 0; i < builtin_getVertices(edges) ; i++) {
      init(i);
    };

    //Build the segmented graph from graph
    int num_places = omp_get_num_places();
    int numSegments = num_places * 5;
    int segmentRange = (edges.num_nodes() + numSegments) / numSegments;
    GraphSegments<int,int>* graphSegments = new GraphSegments<int,int>(numSegments, num_places);
    BuildPullSegmentedGraphs(&edges, graphSegments, segmentRange);

    // Build and initialize per socket local buffer
    localIDs = new int*[num_places];
    for (int socketId = 0; socketId < num_places; socketId++) {
      localIDs[socketId] = (int *)numa_alloc_onnode(sizeof(int) * edges.num_nodes(), socketId);
#pragma omp parallel for
      for (int n = 0; n < edges.num_nodes(); n++) {
	localIDs[socketId][n] = n;
      }
    }

    omp_set_nested(1);

    Timer numa_timer;
    numa_timer.Start();

    while ( (builtin_getVertexSetSize(frontier) ) != ((0) ))
    {
      frontier = edgeset_apply_hybrid_dense_parallel_deduplicatied_from_vertexset_with_frontier(edges, frontier, updateEdge, updateEdge_push_ver, graphSegments); 
    }

#ifdef COUNT
    unordered_map<int, bool> labels = {};
    int cnt = 0;
    for (int i = 0; i < edges.num_nodes();i++){
      if (labels.find(IDs[i]) == labels.end()) {
	cnt++;
	labels[IDs[i]] = true;
      }
    }
    cout << "Number of connected components=" << cnt << endl;
#endif

    numa_timer.Stop();
    PrintTime("NUMA Time", numa_timer.Seconds());
    for (int socketId = 0; socketId < num_places; socketId++)
      numa_free(localIDs[socketId], sizeof(int) * edges.num_nodes());
    delete[] localIDs;

    float elapsed_time = stopTimer() ;
    std::cout << "elapsed time: "<< std::endl;
    std::cout << elapsed_time<< std::endl;
  }
};

