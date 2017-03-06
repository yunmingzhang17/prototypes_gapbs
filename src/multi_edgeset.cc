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


using namespace std;



int main(int argc, char* argv[]) {
    CLApp cli(argc, argv, "multi_edgeset app");
//    if (!cli.ParseArgs())
//        return -1;



    //load in the first graph
    Builder b(cli);
    Graph trust_graph = b.MakeGraph("../test/graphs/filmtrust/trust.el");

    //load in the second graph. Both graphs should have consistent vertex IDs.
    //TODO: fix the weight type from int32_t to float
    WeightedBuilder wb(cli);
    WGraph ratings_graph = wb.MakeGraph("../test/graphs/filmtrust/ratings.wel");

    //pick a series of starting points
    SourcePicker<Graph> sp(trust_graph);


//    auto BFSBound = [&sp] (const Graph &g) { return DOBFS(g, sp.PickNext()); };
//    SourcePicker<Graph> vsp(g, cli.start_vertex());
//    auto VerifierBound = [&vsp] (const Graph &g, const pvector<NodeID> &parent) {
//        return BFSVerifier(g, vsp.PickNext(), parent);
//    };
//    BenchmarkKernel(cli, g, BFSBound, PrintBFSStats, VerifierBound);
    return 0;
}
