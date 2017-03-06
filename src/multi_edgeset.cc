//
// Created by Yunming Zhang on 3/5/17.
//

#include <iostream>
#include <vector>

#include "benchmark.h"
#include "bitmap.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "platform_atomics.h"
#include "pvector.h"
#include "sliding_queue.h"
#include "timer.h"

typedef float WeightFloatT;
typedef NodeWeight<NodeID, WeightFloatT> WFloatNode;
typedef CSRGraph<NodeID, WFloatNode> WFloatGraph;
typedef BuilderBase<NodeID, WFloatNode, WeightFloatT > WeightedFloatBuilder;


using namespace std;


pvector<NodeID> BuildTrustCircle(const Graph &trust_graph, NodeID source){
    pvector<NodeID> trust_circle;
    for (NodeID ngh : trust_graph.out_neigh(source)){
        trust_circle.push_back(ngh);
        cout << "ngh: " << ngh << endl;
    }
    return trust_circle;
}

pvector<NodeID> Recommend(const WFloatGraph &ratings_graph, pvector<NodeID> &trust_circle){
    pvector<NodeID> items;
    for (NodeID influencer : trust_circle){
        for (WFloatNode item : ratings_graph.out_neigh(influencer)){
            if (item.w > 3){
                items.push_back(item.v);
            }
        }
    }
    return items;
}


pvector<NodeID> DoRecommendation(const Graph &trust_graph, const WFloatGraph &ratings_graph, NodeID source){
    PrintStep("Source", static_cast<int64_t>(source));
    Timer t;
    t.Start();
    pvector<NodeID> trust_circle = BuildTrustCircle(trust_graph, source);
    t.Stop();
    PrintStep("Build Circle of Trust", t.Seconds());

    t.Start();
    Recommend(ratings_graph, trust_circle);
    t.Stop();
    PrintStep("Recommend", t.Seconds());

    return trust_circle;
}

int main(int argc, char* argv[]) {
    CLApp cli(argc, argv, "multi_edgeset app");
//    if (!cli.ParseArgs())
//        return -1;



    //load in the first graph
    Builder b(cli);
    Graph trust_graph = b.MakeGraph("../test/graphs/filmtrust/trust.el");

    //load in the second graph. Both graphs should have consistent vertex IDs.
    //TODO: fix the weight type from int32_t to float
    WeightedFloatBuilder wb(cli);
    WFloatGraph ratings_graph = wb.MakeGraph("../test/graphs/filmtrust/ratings.wel");

    //pick a series of starting points
    SourcePicker<Graph> sp(trust_graph);
    DoRecommendation(trust_graph, ratings_graph, sp.PickNext());

//    auto BFSBound = [&sp] (const Graph &g) { return DOBFS(g, sp.PickNext()); };
//    SourcePicker<Graph> vsp(g, cli.start_vertex());
//    auto VerifierBound = [&vsp] (const Graph &g, const pvector<NodeID> &parent) {
//        return BFSVerifier(g, vsp.PickNext(), parent);
//    };
//    BenchmarkKernel(cli, g, BFSBound, PrintBFSStats, VerifierBound);
    return 0;
}
