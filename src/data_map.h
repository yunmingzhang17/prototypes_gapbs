//
// Created by Yunming Zhang on 3/13/17.
//

#ifndef MIGRA_GAPBS_DATA_MAP_H
#define MIGRA_GAPBS_DATA_MAP_H


#include <unordered_map>
#include "pvector.h"

using namespace std;

template <typename T, typename V>
class DataMap {

public:
    DataMap (){}
    virtual ~DataMap(){}
    virtual void init(V init_val) = 0;
    virtual void find_and_add(T index, V val) = 0;
    virtual int num_updated() = 0;
    virtual pvector<T> get_top_keys(int top_count) = 0;
    virtual pvector<T> get_updated_keys() = 0;
};

template <typename T, typename V>
class UnorderedMapDataMap : public DataMap <T,V> {

private:
    V init_val_;
    unordered_map<T, V> data_map_;
    int num_updated_;

public:
    UnorderedMapDataMap<T,V> (size_t init_size) : DataMap<T,V> (){
        num_updated_ = 0;
    }

    virtual ~UnorderedMapDataMap() {
        //delete data_map_;
    }

    virtual void init(V init_val) {
        init_val_ = init_val;
    }

    virtual void find_and_add(T index, V val) {
        // Data Map does not contain the entry
        if (data_map_.find(index) != data_map_.end()) {
            data_map_[index] = init_val_ + val;
            num_updated_++;
        } else {
            data_map_[index] = data_map_[index] + val;
        }
    }

    virtual int num_updated() {
        return num_updated_;
    }

    virtual pvector<T> get_top_keys(int top_count) {
        pvector<T> buffer(num_updated_);
        size_t output_count = top_count < num_updated_ ? top_count : num_updated_;
        pvector<T> output_buffer(output_count);

        for (auto kv : data_map_){
            buffer.push_back(kv.first);
        }
        sort(buffer.begin(), buffer.end(), [this] (T const& a, T const& b) { return data_map_[a] > data_map_[b];});
        for (int i = 0; i < output_count; i++){
            output_buffer.push_back(buffer[i]);
        }
        return output_buffer;
    }

    virtual pvector<T> get_updated_keys() {
        pvector<T> buffer(num_updated_);
        for (auto kv : data_map_){
            buffer.push_back(kv.first);
        }
        return buffer;
    }
};

#endif //MIGRA_GAPBS_DATA_MAP_H
