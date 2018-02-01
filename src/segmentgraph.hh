#include <math.h>
#include <vector>
#include <assert.h>
#include <numa.h>
#include <omp.h>

#include "graph.h"
#include "pvector.h"

using namespace std;

template <class DataT, class Vertex>
struct SegmentedGraph 
{
  int *graphId;
  int *edgeArray;
  int *vertexArray;
  struct SegmentedGraph<DataT, Vertex> *next;
  int numVertices;
  int numEdges;
  bool allocated;

  // debug purpose
  int segmentId;

private:
  int lastLocalIndex;
  Vertex lastVertex;
  Vertex lastEdgeIndex;

public:
  SegmentedGraph(int _id) : segmentId(_id)
  {
    allocated = false;
    next = NULL;
    numVertices = 0;
    numEdges = 0;
    lastVertex = -1;
    lastEdgeIndex = 0;
    lastLocalIndex = 0;
  }

  ~SegmentedGraph()
  {
    numa_free(graphId, sizeof(int) * numVertices);
    numa_free(edgeArray, sizeof(int) * numEdges);
    numa_free(vertexArray, sizeof(int) * (numVertices + 1));
  }


  void allocate(int segment_id)
  {
    int node = segment_id % omp_get_num_places();;

    vertexArray = (int *)numa_alloc_onnode(sizeof(int) * (numVertices + 1), node); // start,end of last
    vertexArray[numVertices] = numEdges;
    edgeArray = (int *)numa_alloc_onnode(sizeof(int) * numEdges, node);
    graphId = (int *)numa_alloc_onnode(sizeof(int) * numVertices, node);
    allocated = true;
  }

  /**
   * Count how many edges we need.
   * @v: dst vertex in pull direction and src vertex in push direction
   **/
  inline
  void countEdge(Vertex v)
  {
    if (v != lastVertex) {
      numVertices++;
      lastVertex = v;
    }
    numEdges++;
  }

  /**
   * Add new edge to each subgraph
   * @v: src in pull direction, dst in push direction
   * @e: dst in pull direction, src in push direction
   **/
  inline
  void addEdge(Vertex toVertexArray, Vertex toEdgeArray)
  {
    if (toVertexArray != lastVertex) {
      // a new vertex going to the same partition                                   
      // must be sorted                                                             
      //assert(e > lastVertex);
      lastVertex = toVertexArray;
      graphId[lastLocalIndex] = toVertexArray;
      vertexArray[lastLocalIndex++] = lastEdgeIndex;
    }                                                                               
    edgeArray[lastEdgeIndex++] = toEdgeArray;
  }

  void print(){
    assert(allocated == true);
    cout << "Segmented Graph numVertices: " << numVertices << " numEdges: " << numEdges  << endl;
    cout << "Vertex Array: " << endl;
    for (int i = 0; i < numVertices; i++){
      cout << " " << vertexArray[i];
    }
    cout << endl;

    cout << "Edge Array: " << endl;
    for (int i = 0; i < numEdges; i++){
      cout << " " << edgeArray[i];
    }
    cout << endl;

    cout << "GraphId Array: " << endl;
    for (int i = 0; i < numVertices; i++){
      cout << " " << graphId[i];
    }
    cout << endl;
  }
};

template <class DataT, class Vertex>
struct GraphSegments 
{
  int numSegments;
  int numGroups;
  SegmentedGraph<DataT,Vertex> **groups;
  SegmentedGraph<DataT,Vertex> **unprocessedGroups;
  vector<SegmentedGraph<DataT,Vertex>*> segments;
  
  GraphSegments(int _numSegments, int _numGroups): numSegments(_numSegments), numGroups(_numGroups)
  {
    groups = new SegmentedGraph<DataT,Vertex>*[numGroups]();
    unprocessedGroups = new SegmentedGraph<DataT,Vertex>*[numGroups]();

    for (int i = 0; i < numSegments; i++){
      auto sg = new SegmentedGraph<DataT, Vertex>(i);
      segments.push_back(sg);
      int groupId = i % numGroups;
      sg->next = groups[groupId];
      groups[groupId] = sg;
    }    
  }

  ~GraphSegments(){
    //deallocate every graph segment
   for (int i=0; i<numSegments; i++){
     delete segments[i];
   }
   delete[] groups;
   delete[] unprocessedGroups;
  }

  void resetGroups() {
    // In every iteration, unprocessedGroups should be set to groups
    memcpy((void *)unprocessedGroups, (void *)groups, sizeof(*groups) * numGroups);
  }

  void allocate() {
    for (int i = 0; i<numSegments; i++){
      segments[i]->allocate(i);
    }
  }

  SegmentedGraph<DataT, Vertex> * getSegmentedGraph(int id){
    return segments[id];
  }

  SegmentedGraph<DataT, Vertex> *getSegmentedGraphFromGroup(int groupId) {

    // Try groupId's group first. If empty, try other grouops
    auto sg = getHeadOrNull(groupId);
    if (sg)
      return sg;
    int victimId = (groupId + 1) % numGroups;
    while (victimId != groupId) {
      sg = getHeadOrNull(victimId);
      if (sg)
	return sg;
      victimId = (victimId + 1) % numGroups;
    }
    return NULL;
  }

  SegmentedGraph<DataT, Vertex> *getHeadOrNull(int groupId) {
    auto head = unprocessedGroups[groupId];
    while (head) {
      auto next = head->next;
      if (__sync_bool_compare_and_swap(&unprocessedGroups[groupId], head, next))
	return head;
      head = unprocessedGroups[groupId];
    }
    return NULL;
  }
};

template <class DataT, class Vertex>
void BuildPullSegmentedGraphs(const Graph* originalGraph, GraphSegments<DataT,Vertex> * graphSegments, int segmentRange)
{
  //Go through the original graph and count the number of target vertices and edges for each segment
  for (NodeID v : originalGraph->vertices()){
    for (NodeID u : originalGraph->in_neigh(v)){
      int segment_id = u/segmentRange;
      graphSegments->getSegmentedGraph(segment_id)->countEdge(v);
    }
  }

  //Allocate each segment
  graphSegments->allocate();

  //Add the edges for each segment
  for (NodeID v : originalGraph->vertices()){
    for (NodeID u : originalGraph->in_neigh(v)){
      int segment_id = u/segmentRange;
      graphSegments->getSegmentedGraph(segment_id)->addEdge(v, u);
    }
  }
}

template <class DataT, class Vertex>
void BuildPushSegmentedGraphs(const Graph* originalGraph, GraphSegments<DataT,Vertex> * graphSegments, int segmentRange)
{
  //Go through the original graph and count the number of target vertices and edges for each segment
  for (NodeID u : originalGraph->vertices()){
    for (NodeID v : originalGraph->out_neigh(u)){
      int segment_id = v/segmentRange;
      graphSegments->getSegmentedGraph(segment_id)->countEdge(u);
    }
  }

  //Allocate each segment
  graphSegments->allocate();

  //Add the edges for each segment
  for (NodeID u : originalGraph->vertices()){
    for (NodeID v : originalGraph->out_neigh(u)){
      int segment_id = v/segmentRange;
      graphSegments->getSegmentedGraph(segment_id)->addEdge(u, v);
    }
  }
}

