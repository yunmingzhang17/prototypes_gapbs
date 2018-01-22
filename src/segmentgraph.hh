#include <math.h>
#include <vector>
#include <assert.h>
#include <numa.h>

#include "graph.h"
#include "pvector.h"

using namespace std;

static int num_numa_node = numa_num_configured_nodes();

template <class DataT, class Vertex>
struct SegmentedGraph 
{
  float *incoming_total;
  int *graphId;
  int *edgeArray;
  int *dstVertexArray;
  int numDstVertices;
  int numEdges;
  int numVertices;
  bool allocated;

private:
  int lastLocalIndex;
  Vertex lastVertex;
  Vertex lastEdgeIndex;

public:
  SegmentedGraph(int _numVertices) : numVertices(_numVertices)
  {
    allocated = false;
    numDstVertices = 0;
    numEdges = 0;
    lastVertex = -1;
    lastEdgeIndex = 0;
    lastLocalIndex = 0;
  }

  ~SegmentedGraph()
  {
    numa_free(incoming_total, sizeof(float) * numVertices);
    numa_free(graphId, sizeof(int) * numDstVertices);
    numa_free(edgeArray, sizeof(int) * numEdges);
    numa_free(dstVertexArray, sizeof(int) * numDstVertices);
    //delete[] graphId;
    //delete[] edgeArray;
    //delete[] dstVertexArray;
    //delete[] incoming_total;
  }


  void allocate(int segment_id)
  {
    int node = segment_id % num_numa_node;

    incoming_total = (float *)numa_alloc_onnode(sizeof(float) * numVertices, node);
    //incoming_total = new float[numVertices];

    #pragma omp parallel for
    for (int i = 0; i < numVertices; i++)
      incoming_total[i] = 0;

    dstVertexArray = (int *)numa_alloc_onnode(sizeof(int) * (numDstVertices + 1), node); // start,end of last
    //dstVertexArray = new int[numDstVertices + 1];
    dstVertexArray[numDstVertices] = numEdges;
    edgeArray = (int *)numa_alloc_onnode(sizeof(int) * numEdges, node);
    //edgeArray = new int[numEdges];
    graphId = (int *)numa_alloc_onnode(sizeof(int) * numDstVertices, node);
    //graphId = new int[numVertices];
    allocated = true;
  }

  //countIncomingEdge counts how many edges we need and create a mapping to global id
  inline
  void countIncomingEdge(Vertex src, Vertex dst)
  {
    bool newVertex = false;
    if (dst != lastVertex) {
      numDstVertices++;
      lastVertex = dst;
      newVertex = true;
    }
    numEdges++;
  }

  //addIncomingEdge adds new edge to each subgraph, pulling from src to dst
  inline
  void addIncomingEdge(Vertex src, Vertex dst)
  {
    if (dst != lastVertex) {
      // a new vertex going to the same partition                                   
      // must be sorted                                                             
      //assert(dst > lastVertex);
      lastVertex = dst;
      graphId[lastLocalIndex] = dst;
      dstVertexArray[lastLocalIndex++] = lastEdgeIndex;
    }                                                                               
    edgeArray[lastEdgeIndex++] = src;
  }

  void print(){
    assert(allocated == true);
    cout << "Segmented Graph numDstVertices: " << numDstVertices << " numEdges: " << numEdges  << endl;
    cout << "DstVertex Array: " << endl;
    for (int i = 0; i < numDstVertices; i++){
      cout << " " << dstVertexArray[i];
    }
    cout << endl;

    cout << "Edge Array: " << endl;
    for (int i = 0; i < numEdges; i++){
      cout << " " << edgeArray[i];
    }
    cout << endl;

    cout << "GraphId Array: " << endl;
    for (int i = 0; i < numDstVertices; i++){
      cout << " " << graphId[i];
    }
    cout << endl;
  }
};

template <class DataT, class Vertex>
struct GraphSegments 
{
  int numSegments;
  vector<SegmentedGraph<DataT,Vertex>*> segments;
  
  GraphSegments(int _numSegments, int numVertices): numSegments(_numSegments)
  {
    //alocate each graph segment
    for (int i=0; i<numSegments; i++){
      segments.push_back(new SegmentedGraph<int, int>(numVertices));
    }
  }

  ~GraphSegments(){
    //deallocate every graph segment
   for (int i=0; i<numSegments; i++){
     delete segments[i];
   }
  }

  void allocate() {
    for (int i = 0; i<numSegments; i++){
      segments[i]->allocate(i);
    }
  }

  SegmentedGraph<DataT, Vertex> * getSegmentedGraph(int id){
    return segments[id];
  }
};

/**
 * Build Graph Segments with input Graph, output graphSegments and specified number of segments. This has to work with GAPBS graphs now
 */
template <class DataT, class Vertex>
void BuildCacheSegmentedGraphs(const Graph* originalGraph, GraphSegments<DataT,Vertex> * graphSegments, int segmentRange)
{
  //Go through the original graph and count the number of target vertices and edges for each segment
  for (NodeID v : originalGraph->vertices()){
    for (NodeID u : originalGraph->in_neigh(v)){
      int segment_id = u/segmentRange;
      graphSegments->getSegmentedGraph(segment_id)->countIncomingEdge(u, v);
    }
  }

  //Allocate each segment
  graphSegments->allocate();

  //Add the edges for each segment
  for (NodeID v : originalGraph->vertices()){
    for (NodeID u : originalGraph->in_neigh(v)){
      int segment_id = u/segmentRange;
      graphSegments->getSegmentedGraph(segment_id)->addIncomingEdge(u, v);
    }
  }
}
