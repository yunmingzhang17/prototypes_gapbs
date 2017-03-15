//
// Created by Yunming Zhang on 3/5/17.
//

#include <iostream>
#include <vector>
#include <unordered_map>

#include "benchmark.h"
#include "bitmap.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "platform_atomics.h"
#include "pvector.h"
#include <queue>
#include <curses.h>
#include "timer.h"
#include "sliding_queue.h"
#include "data_map.h"


#define DEBUG_DETAILS

#define USE_STL_HASHMAP
#define USE_GOOGLE_HASHMAP
#define USE_ARRAY

//typedef float WeightFloatT;
//typedef NodeWeight<NodeID, WeightFloatT> WFloatNode;
//typedef CSRGraph<NodeID, WFloatNode> WFloatGraph;
//typedef BuilderBase<NodeID, WFloatNode, WeightFloatT > WeightedFloatBuilder;
//

using namespace std;


vector<NodeID> BuildTrustCircle(const Graph &trust_graph, NodeID source){
    vector<NodeID> trust_circle;
    Bitmap visited(trust_graph.num_nodes());
    visited.reset();
    vector<NodeID>* to_visit = new vector<NodeID>;
    to_visit->reserve(trust_graph.num_nodes());
    vector<NodeID>* to_visit_next = new vector<NodeID>;
    to_visit_next->reserve(trust_graph.num_nodes());

    to_visit->push_back(source);


    //three hops of all the unique nodes
    int max_hop = 3;
    int cur_hop = 0;

    while (!to_visit->empty() && cur_hop < max_hop){
#ifdef DEBUG_DETAILS
        cout << "queue size: " << to_visit->size() << endl;
#endif
        for (auto active_vertex : *to_visit){
            for (NodeID ngh : trust_graph.out_neigh(active_vertex)){
#ifdef DEBUG_DETAILS
                //cout << "active_vertex: " << active_vertex << " ngh: " << ngh << endl;
#endif
                if (!visited.get_bit(ngh)) {
                    visited.set_bit(ngh);
                    to_visit_next->push_back(ngh);
                    trust_circle.push_back(ngh);
                }
            }
        }
        cur_hop++;
        to_visit = to_visit_next;
        to_visit_next = new vector<NodeID>();
    }

    return trust_circle;
}

pvector<NodeID> RecommendDataMap(const WGraph &ratings_graph, vector<NodeID> &trust_circle, DataMap<NodeID, int> &data_map){

    data_map.init(0);

    for (NodeID influencer : trust_circle){
        //A node that is present in the social graph, but is larger than the largest people node in ratings graph
        if (influencer >= ratings_graph.num_nodes()){
            continue;
        }
        for (WNode item : ratings_graph.out_neigh(influencer)){
            if (item.w > 3){
                data_map.find_and_add(item.v, 1);
            }
        }
    }

#ifdef DEBUG_DETAILS
    cout << "count map size: " << data_map.num_updated() << endl;
    cout << "top counts: " << endl;
    pvector<NodeID> items = data_map.get_top_keys(5);
    for (int i = 0; i < items.size(); i++){
        cout << "item: " << items[i] << " count: " << data_map.find(items[i]) << endl;
    }
#endif

    return data_map.get_updated_keys();

}

pvector<NodeID> RecommendHashMap(const WGraph &ratings_graph, vector<NodeID> &trust_circle){
    pvector<NodeID> items;
    unordered_map<NodeID,int32_t> count_map;

    for (NodeID influencer : trust_circle){
        //A node that is present in the social graph, but is larger than the largest people node in ratings graph
        if (influencer >= ratings_graph.num_nodes()){
            continue;
        }
        for (WNode item : ratings_graph.out_neigh(influencer)){
            if (item.w > 3){
                //items.push_back(item.v);
                if (count_map.find(item.v) != count_map.end()){
                  count_map[item.v] = count_map[item.v]++;
                }else{
                    count_map[item.v] = 1;
                }
            }
        }
    }

    items.reserve(count_map.size());
    for (auto kv : count_map){
        items.push_back(kv.first);
    }
    sort(items.begin(), items.end(), [&count_map] (NodeID const& a, NodeID const& b) {
        bool output = FALSE;
        if (count_map[a] > count_map[b]){
            output = TRUE;
        } else if (count_map[a] == count_map[b]){
            if (a > b){
                output = TRUE;
            }
        }
        return output;
    });

#ifdef DEBUG_DETAILS
    cout << "count map size: " << count_map.size() << endl;
    cout << "top counts: " << endl;
    int count = 5 < count_map.size() ? 5 : count_map.size();
    for (int i = 0; i < count; i++){
        cout << "item: " << items[i] << " count: " << count_map[items[i]] << endl;
    }
#endif

    return items;
}

pvector<NodeID> RecommendArray(const WGraph &ratings_graph, vector<NodeID> &trust_circle, int32_t num_items){
    pvector<NodeID> items;
    vector<int> count_map(num_items);
    for (int i = 0; i < num_items; i++){
        count_map[i] = 0;
    }

    for (NodeID influencer : trust_circle){
        //A node that is present in the social graph, but is larger than the largest people node in ratings graph
        if (influencer >= ratings_graph.num_nodes()){
            continue;
        }

        for (WNode item : ratings_graph.out_neigh(influencer)){
            if (item.w > 3){
                //items.push_back(item.v);
                count_map[item.v]++;
            }
        }
    }

    //items.reserve(count_map.size());
    for (int i = 0; i < num_items; i++){
        if (count_map[i] > 0) {
            items.push_back(i);
        }
    }
    sort(items.begin(), items.end(), [&count_map] (NodeID const& a, NodeID const& b) {
        bool output = FALSE;
        if (count_map[a] > count_map[b]){
            output = TRUE;
        } else if (count_map[a] == count_map[b]){
            if (a > b){
                output = TRUE;
            }
        }
        return output;
    });

#ifdef DEBUG_DETAILS
    cout << "count map size: " << count_map.size() << endl;
    cout << "top counts: " << endl;
    int count = 5 < items.size() ? 5 : items.size();
    for (int i = 0; i < count; i++){
        cout << "item: " << items[i] << " count: " << count_map[items[i]] << endl;
    }
#endif

    return items;
}

vector<NodeID> ParBuildTrustCircle(const Graph &trust_graph, NodeID source){
    vector<NodeID> trust_circle;
    //Bitmap visited(trust_graph.num_nodes());
    pvector<bool> visited(trust_graph.num_nodes());
    SlidingQueue<NodeID> queue(trust_graph.num_nodes());
    queue.push_back(source);
    queue.slide_window();

    //three hops of all the unique nodes
    int max_hop = 3;
    int cur_hop = 0;


    for(int i = 0; i < trust_graph.num_nodes(); i++){
        visited[i] = FALSE;
    }

    while (!queue.empty() && cur_hop < max_hop){
#ifdef DEBUG_DETAILS
        cout << "queue size: " << queue.size() << endl;
#endif
#pragma omp parallel
        {
            QueueBuffer<NodeID> lqueue(queue);
#pragma omp for
            for (auto q_iter = queue.begin(); q_iter < queue.end(); q_iter++) {
                NodeID u = *q_iter;
                for (NodeID ngh : trust_graph.out_neigh(u)) {
#ifdef DEBUG_DETAILS
                    //cout << "active_vertex: " << active_vertex << " ngh: " << ngh << endl;
#endif
                    if (!visited[ngh]) {
                        if (compare_and_swap(visited[ngh], FALSE, TRUE)){
                            lqueue.push_back(ngh);
                        }
                    }
                }
            }
            lqueue.flush();
        }
        queue.slide_window();

        for (auto node : queue){
            trust_circle.push_back(node);
        }

        cur_hop++;
    }
    cout << "trust circle size: " << trust_circle.size() << endl;
    return trust_circle;
}

// A version that uses local unordered_maps in each thread and later merge them
pvector<NodeID> ParRecommendLocalMap(const WGraph &ratings_graph, vector<NodeID> &trust_circle){
    pvector<NodeID> items;
    unordered_map<NodeID,int32_t> count_map;
    int top_count = 10;

    #pragma omp parallel
    {
        unordered_map<NodeID,int32_t> local_count_map;
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < trust_circle.size(); i++) {
            NodeID influencer = trust_circle[i];
            //A node that is present in the social graph, but is larger than the largest people node in ratings graph
            if (influencer >= ratings_graph.num_nodes()){
                continue;
            }
            for (WNode item : ratings_graph.out_neigh(influencer)) {
                if (item.w > 3) {
                    //items.push_back(item.v);
                    if (local_count_map.find(item.v) != local_count_map.end()) {
                        local_count_map[item.v] = local_count_map[item.v]++;
                    } else {
                        local_count_map[item.v] = 1;
                    }
                }
            }
        }

        #pragma omp critical
        {
            for (auto kv : local_count_map) {
                if (count_map.find(kv.first) != count_map.end())
                    count_map[kv.first] += local_count_map[kv.first];
                else
                    count_map[kv.first] = local_count_map[kv.first];
            }
        }

    }


    items.reserve(count_map.size());
    for (auto kv : count_map){
        items.push_back(kv.first);
    }
    sort(items.begin(), items.end(), [&count_map] (NodeID const& a, NodeID const& b) {
        bool output = FALSE;
        if (count_map[a] > count_map[b]){
            output = TRUE;
        } else if (count_map[a] == count_map[b]){
            if (a > b){
                output = TRUE;
            }
        }
        return output;
    });

#ifdef DEBUG_DETAILS
    cout << "count map size: " << count_map.size() << endl;
    cout << "top counts: " << endl;
    int count = 5 < count_map.size() ? 5 : count_map.size();
    for (int i = 0; i < count; i++){
        cout << "item: " << items[i] << " count: " << count_map[items[i]] << endl;
    }
#endif

    return items;
}


// A version that uses hash locked hashmap
pvector<NodeID> ParRecommendArray(const WGraph &ratings_graph, vector<NodeID> &trust_circle, int32_t num_items){
    pvector<NodeID> items;
    vector<int> count_map(num_items);

    #pragma omp parallel for
    for (int i = 0; i < num_items; i++){
        count_map[i] = 0;
    }

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < trust_circle.size(); i++){
        NodeID influencer = trust_circle[i];
        //A node that is present in the social graph, but is larger than the largest people node in ratings graph
        if (influencer >= ratings_graph.num_nodes()){
            continue;
        }
        for (WNode item : ratings_graph.out_neigh(influencer)){
            if (item.w > 3){
                //items.push_back(item.v);
                fetch_and_add(count_map[item.v],1);
            }
        }
    }

    //items.reserve(count_map.size());
    for (int i = 0; i < num_items; i++){
        if (count_map[i] > 0) {
            items.push_back(i);
        }
    }
    sort(items.begin(), items.end(), [&count_map] (NodeID const& a, NodeID const& b) {
        bool output = FALSE;
        if (count_map[a] > count_map[b]){
            output = TRUE;
        } else if (count_map[a] == count_map[b]){
            if (a > b){
                output = TRUE;
            }
        }
        return output;
    });

#ifdef DEBUG_DETAILS
    cout << "count map size: " << count_map.size() << endl;
    cout << "top counts: " << endl;
    int count = 5 < items.size() ? 5 : items.size();
    for (int i = 0; i < count; i++){
        cout << "item: " << items[i] << " count: " << count_map[items[i]] << endl;
    }
#endif

    return items;
}


vector<NodeID> DoRecommendation(const Graph &trust_graph, const WGraph &ratings_graph, NodeID source, int num_items, bool is_parallel){
    PrintStep("Source", static_cast<int64_t>(source));
    vector<NodeID> top_items(5);
    pvector<NodeID> items;

    if (is_parallel){
        //parallel execution
        Timer t;
        t.Start();
        vector<NodeID> trust_circle = ParBuildTrustCircle(trust_graph, source);
        t.Stop();
        PrintStep("Parallel Build Circle of Trust", t.Seconds());

#ifdef DEBUG_DETAILS
        //    for (auto trustee : trust_circle){
        //        cout << "trustee: " << trustee << endl;
        //    }
            cout << "number of trustees: " << trust_circle.size() << endl;
#endif

#ifdef USE_STL_HASHMAP
        t.Start();
        items = ParRecommendLocalMap(ratings_graph, trust_circle);
        t.Stop();
        PrintStep("Parallel Hash Map Recommendation", t.Seconds());
#endif

#ifdef USE_ARRAY
        t.Start();
        items = ParRecommendArray(ratings_graph, trust_circle, num_items);
        t.Stop();
        PrintStep("Parallel Array Recommendation", t.Seconds());
#endif

    } else {
        //serial execution
        Timer t;
        t.Start();
        vector<NodeID> trust_circle = BuildTrustCircle(trust_graph, source);
        t.Stop();
        PrintStep("Build Circle of Trust", t.Seconds());

#ifdef DEBUG_DETAILS
        //    for (auto trustee : trust_circle){
        //        cout << "trustee: " << trustee << endl;
        //    }
            cout << "number of trustees: " << trust_circle.size() << endl;
#endif

#ifdef USE_STL_HASHMAP
        t.Start();
        //items = RecommendHashMap(ratings_graph, trust_circle);
        UnorderedMapDataMap<NodeID, int> unordered_data_map(num_items);
        items = RecommendDataMap(ratings_graph, trust_circle, unordered_data_map);
        t.Stop();
        PrintStep("STL Serial Hash Map based Recommendation", t.Seconds());
#endif

#ifdef USE_GOOGLE_HASHMAP
        t.Start();
        GoogleHashDataMap<NodeID, int> google_data_map(num_items);
        items = RecommendDataMap(ratings_graph, trust_circle, google_data_map);
        t.Stop();
        PrintStep("Google Serial Hash Map based Recommendation", t.Seconds());
#endif

#ifdef USE_ARRAY
        t.Start();
        ArrayDataMap<NodeID, int> array_data_map(num_items);
        items = RecommendDataMap(ratings_graph, trust_circle, array_data_map);
        //items = RecommendArray(ratings_graph, trust_circle, num_items);
        t.Stop();
        PrintStep("Serial Array based Recommendation", t.Seconds());
#endif

    }

    int count = 5 < items.size() ? 5 : items.size();
    for (int i = 0; i < count; i++){
        top_items[i] = items[i];
    }
    return top_items;
}


void PrintStats(const Graph &g, vector<NodeID> &items) {
    for (int i = 0; i < items.size(); i++){
        cout << "item: " << items[i]  << endl;
    }
}

int main(int argc, char* argv[]) {
    CLMultiEdgeSet cli(argc, argv, "multi_edgeset app");
    if (!cli.ParseArgs())
        return -1;



    //load in the first graph
    Builder b(cli);
    //Graph trust_graph = b.MakeGraph("../test/graphs/filmtrust/trust_revised.el");
    Graph trust_graph = b.MakeGraph();

    //load in the second graph. Both graphs should have consistent vertex IDs.
    WeightedBuilder wb(cli);
    //WFloatGraph ratings_graph = wb.MakeGraph("../test/graphs/filmtrust/ratings.wel");
    WGraph ratings_graph = wb.MakeSecondGraph();

    //pick a series of starting points
    SourcePicker<Graph> sp(trust_graph, cli.start_vertex());
//    for (int i = 0; i < 10; i++){
//        NodeID start = sp.PickNext();
//        DoRecommendation(trust_graph, ratings_graph, start, cli.num_items(), FALSE);
//        DoRecommendation(trust_graph, ratings_graph, start, cli.num_items(), TRUE);
//    }

    int num_items = cli.num_items() + 1;//deal with potential ID issue

    auto SimpleRecBound = [&sp, &ratings_graph, &num_items] (const Graph &g) {
        return DoRecommendation(g, ratings_graph, sp.PickNext(), num_items, TRUE); };
    SourcePicker<Graph> vsp(trust_graph, cli.start_vertex());

    auto VerifierBound = [&vsp, &ratings_graph, &num_items] (const Graph &g, const vector<NodeID> &parallel_top_items) {
        vector<NodeID> serial_top_items = DoRecommendation(g, ratings_graph, vsp.PickNext(), num_items, FALSE);
        bool output = TRUE;

        if (serial_top_items.size() != parallel_top_items.size()){
            output = FALSE;
            cout << " size of output do not match " << endl;
        }


        for (int i = 0; i < serial_top_items.size(); i++){
            if (serial_top_items[i] != parallel_top_items[i]){
                output = FALSE;
                cout << " items do not match " << endl;
                cout << "serial: " << serial_top_items[i] << endl;
                cout << "parallel: " << parallel_top_items[i] << endl;
                break;
            }
        }

        return output;
    };
    BenchmarkKernel(cli, trust_graph, SimpleRecBound, PrintStats, VerifierBound);
    return 0;
}
