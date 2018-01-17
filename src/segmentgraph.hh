#include <math.h>
#include <vector>
#include <assert.h>
#include <numa.h>

#include "graph.h"
#include "pvector.h"

using namespace std;

int num_numa_node = numa_num_configured_nodes();

template <class DataT, class Vertex>
struct SegmentedGraph 
{
  //pvector<ScoreT> scores;
  float * scores;
  float * outgoing_contrib;
  int totalVCnt;
  int * graphId;
  int * edgeArray;
  int * vertexArray;
  int numVertices;
  int numEdges;
  int numSlices;
  int numElementsPerSlice;
  bool allocated;
  Timer cumulativeTimer;
  DataT * intermediateBuffer;
  DataT ** sliceIdx;

private:
  int lastLocalIndex;
  Vertex lastVertex;
  Vertex lastEdge;

public:
  SegmentedGraph(int _numSlices, int _numElementsPerSlice) : numSlices(_numSlices) , numElementsPerSlice(_numElementsPerSlice)
  {
    allocated = false;
    numVertices = 0;
    numEdges = 0;
    lastVertex = 0;
    lastEdge = 0;
    lastLocalIndex = 0;
  }

  ~SegmentedGraph()
  {
    numa_free(scores, sizeof(float) * totalVCnt);
    numa_free(outgoing_contrib, sizeof(float) * totalVCnt);
    numa_free(graphId, sizeof(int) * numVertices);
    numa_free(edgeArray, sizeof(int) * numVertices);
    numa_free(vertexArray, sizeof(int) * numVertices);
  }


  void allocate(int total_v_cnt, int segment_id)
  {
    float init_val = 1.0f / numVertices;
    int node = segment_id % num_numa_node;

    scores = (float *)numa_alloc_onnode(sizeof(float) * total_v_cnt, node);

    #pragma omp parallel for
    for (int i = 0; i < total_v_cnt; i++)
      scores[i] = init_val;

    outgoing_contrib = (float *)numa_alloc_onnode(sizeof(float) * total_v_cnt, node);

    vertexArray = (int *)numa_alloc_onnode(sizeof(int) * (numVertices + 1), node); // start,end of last
    vertexArray[numVertices] = numEdges;
    edgeArray = (int *)numa_alloc_onnode(sizeof(int) * numEdges, node);
    graphId = (int *)numa_alloc_onnode(sizeof(int) * numVertices, node);
    lastVertex = 0;
    lastEdge = 0;
    lastLocalIndex = 0;
    totalVCnt = total_v_cnt;
    allocated = true;
  }

  //countIncomingEdge counts how many edges we need and create a mapping to global id
  inline
  bool countIncomingEdge(Vertex src, Vertex dst)
  {
    bool newVertex = false;
    if (dst != lastVertex) {
      numVertices++;
      lastVertex = dst;
      newVertex = true;
    }

    numEdges++;
    return newVertex;
  }

  //addIncomingEdge adds new edge to each subgraph, pulling from src to dst
  inline
  bool addIncomingEdge(Vertex src, Vertex dst)
  {
    if (dst != lastVertex) {
      // a new vertex going to the same partition                                   
      // must be sorted                                                             
      //assert(dst > lastVertex);
      lastVertex = dst;
      graphId[lastLocalIndex] = dst;
      vertexArray[lastLocalIndex++] = lastEdge;
    }                                                                               
    edgeArray[lastEdge++] = src;
    return false;
  }

  void print(){

    assert(allocated == true);
    cout << "Blocked Graph numVertices: " << numVertices << " numEdges: " << numEdges  << endl;
    cout << "Vertex Array: " << endl;
    for (int i=0; i<numVertices; i++){
      cout << " " << vertexArray[i];
    }
    cout << endl;

    cout << "GraphId Array: " << endl;
    for (int i=0; i<numVertices; i++){
      cout << " " << graphId[i];
    }
    cout << endl;

    cout << "Edge Array: " << endl;
    for (int i=0; i<numEdges; i++){
      cout << " " << edgeArray[i];
    }
    cout << endl;
  }

};

template <class DataT, class Vertex>
struct GraphSegments 
{
  int numSegments;
  int numSlices;
  int numElementsPerSlice;
  vector<SegmentedGraph<DataT,Vertex>*> segments;
  
  GraphSegments(int _numSegments, int _numSlices, int _numElementsPerSlice): numSegments(_numSegments), numSlices(_numSlices), numElementsPerSlice(_numElementsPerSlice)
  {
    //alocate each graph segment
    for (int i=0; i<numSegments; i++){
      segments.push_back(new SegmentedGraph<int, int>(numSlices, numElementsPerSlice));
    }
  }

  ~GraphSegments(){
    //deallocate every graph segment
   for (int i=0; i<numSegments; i++){
     delete segments[i];
   }
   //delete blocks;
  }

  void allocate(int total_v_cnt) {
    for (int i=0; i<numSegments; i++){
      segments[i]->allocate(total_v_cnt, i);
    }
  }

  SegmentedGraph<DataT, Vertex> * getSegmentedGraph(int id){
    return segments[id];
  }

};




/**
 * Build Graph Segments with input Graph, output graphSegments and  specified number of segments. This has to work with GAPBS graphs now
 */


template <class DataT, class Vertex>
void BuildCacheSegmentedGraphs(const Graph* originalGraph, GraphSegments<DataT,Vertex> * graphSegments,  int numSlices, int numElementsPerSegment, int numElementsPerSlice)
{
  int64_t numVertices = originalGraph->num_nodes();
  int64_t numEdges = originalGraph->num_edges();
  uint64_t * offsets = new uint64_t[numVertices];
  //pvector<SGOffset> offsets = originalGraph->VertexOffsets();
  //int* vertexArray = originalGraph->vertexArray;
  //int* edgeArray = originalGraph->edgeArray;

  //Right now, this is targeting double vertex value (PageRank)
  //int numBlocksPerRow = ((int) sqrt(numBlocks));
  //int blockDim = ceil((double)numVertices/numBlocksPerRow);

  uint64_t curOffset = 0;
  for (NodeID u : originalGraph->vertices()){
    assert(curOffset <= 2*originalGraph->num_edges());
    offsets[u] = curOffset;
    curOffset += originalGraph->out_degree(u);
  }

  //Go through the original graph and count the number of vertices and edges for each block

  for (NodeID u : originalGraph->vertices()){
    //int start = vertexArray[i];
    //int end = (i == numVertices-1? numEdges:vertexArray[i+1]);
    for (NodeID v : originalGraph->out_neigh(u)){
      int segment_id = offsets[v]/numElementsPerSegment;
      //#ifdef DEBUG
      //cout << "edge from u: " << u << " to v: " << v << endl;
      //#endif
      graphSegments->getSegmentedGraph(segment_id)->countIncomingEdge(v,u);
    }
  }

#ifdef DEBUG
  cout << endl;
#endif
  cout << "allocating segments" << endl;
  //Allocate eac block
  graphSegments->allocate(numVertices);

  //Add the edges for each block
  for (NodeID u : originalGraph->vertices()){
    //int start = vertexArray[i];
    //int end = (i == numVertices-1? numEdges:vertexArray[i+1]);              
    for (NodeID v : originalGraph->out_neigh(u)){
      //int ngh = edgeArray[j];                                                 
      //int blk_id = i/blockDim * numBlocksPerRow + ngh/blockDim;               
      int segment_id = offsets[v]/numElementsPerSegment;
#ifdef DEBUG3
      cout << "dst: " << u << " src: " << v << " blk_id: " << segment_id << endl;
#endif
      graphSegments->getSegmentedGraph(segment_id)->addIncomingEdge(v,u);
    }
  }
}




//Reset the timers in each blocked graph 
template <class DataT, class Vertex>
void ResetBlockTimes(GraphSegments<DataT, Vertex>* graphSegments)
{
  for (auto sg : graphSegments->segments){
    sg->cumulativeTimer.reset();
  }
}

template <class DataT, class Vertex>
void ShowBlockTimes(GraphSegments<DataT, Vertex>* graphSegments) {
  for (int i=0; i<graphSegments->numSegments; i++) {
    auto sg = graphSegments->getSegmentedGraph(i);
    cout << "Segment: " << i
	 << " numVertices: " << sg->numVertices
	 << " numEdges: " << sg->numEdges
	 << " ) took: ";
      sg->cumulativeTimer.printTimes(cout) << endl;
  }
}


