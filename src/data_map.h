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
    virtual V find(T index) = 0;
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

    virtual V find(T index) {
        return data_map_[index];
    }


    virtual void find_and_add(T index, V val) {
        // Data Map does not contain the entry
        if (data_map_.find(index) != data_map_.end()) {
            data_map_[index] = data_map_[index] + val;
        } else {
            data_map_[index] = init_val_ + val;
            num_updated_++;
        }
    }

    virtual int num_updated() {
        return num_updated_;
    }

    virtual pvector<T> get_top_keys(int top_count) {
        pvector<T> buffer(num_updated_);
        size_t output_count = top_count < num_updated_ ? top_count : num_updated_;
        pvector<T> output_buffer(output_count);

        int i = 0;
        for (auto kv : data_map_){
            buffer[i] = kv.first;
            i++;
        }

        sort(buffer.begin(), buffer.end(), [this] (NodeID const& a, NodeID const& b) {
            bool output = FALSE;
            if (data_map_[a] > data_map_[b]){
                output = TRUE;
            } else if (data_map_[a] == data_map_[b]){
                if (a > b){
                    output = TRUE;
                }
            }
            return output;
        });

        for (int i = 0; i < output_count; i++){
            output_buffer[i] = buffer[i];
        }
        return output_buffer;
    }

    virtual pvector<T> get_updated_keys() {
        pvector<T> buffer(num_updated_);
        int i = 0;
        for (auto kv : data_map_){
            buffer[i] = kv.first;
            i++;
        }
        sort(buffer.begin(), buffer.end(), [this] (NodeID const& a, NodeID const& b) {
            bool output = FALSE;
            if (data_map_[a] > data_map_[b]){
                output = TRUE;
            } else if (data_map_[a] == data_map_[b]){
                if (a > b){
                    output = TRUE;
                }
            }
            return output;
        });

        return buffer;
    }
};


template <typename T, typename V>
class ArrayDataMap : public DataMap <T,V> {

private:
    V init_val_;
    pvector<V> data_array_;

public:
    ArrayDataMap<T,V> (size_t init_size) : DataMap<T,V> (){
        data_array_.resize(init_size);
    }

    virtual ~ArrayDataMap() {
        //delete data_map_;
    }

    virtual void init(V init_val) {
        for (int i = 0; i < data_array_.size(); i++){
            data_array_[i] = init_val;
        }
        init_val_ = init_val;
    }

    virtual V find(T index) {
        return data_array_[index];
    }


    virtual void find_and_add(T index, V val) {
        data_array_[index] = data_array_[index] + val;
    }

    virtual int num_updated() {
        int num_updated = 0;
        for (auto val : data_array_){
            if (val != init_val_)
                num_updated++;
        }
        return num_updated;

    }

    virtual pvector<T> get_top_keys(int top_count) {
        pvector<T> buffer(data_array_.size());
        int num_updated = this->num_updated();
        if ( num_updated == 0){
            buffer.resize(0);
            return buffer;
        }

        size_t output_count = top_count < num_updated ? top_count : num_updated;
        pvector<T> output_buffer(output_count);

        for (int i = 0; i < data_array_.size(); i++){
            buffer[i] = i;
        }

        sort(buffer.begin(), buffer.end(), [this] (NodeID const& a, NodeID const& b) {
            bool output = FALSE;
            if (data_array_[a] > data_array_[b]){
                output = TRUE;
            } else if (data_array_[a] == data_array_[b]){
                if (a > b){
                    output = TRUE;
                }
            }
            return output;
        });

        for (int i = 0; i < output_count; i++){
            output_buffer[i] = buffer[i];
        }
        return output_buffer;
    }

    virtual pvector<T> get_updated_keys() {
        pvector<T> buffer(data_array_.size());
        if (this->num_updated() == 0){
            buffer.resize(0);
            return buffer;
        }

        for (int i = 0; i < data_array_.size(); i++){
            buffer[i] = i;
        }

        sort(buffer.begin(), buffer.end(), [this] (NodeID const& a, NodeID const& b) {
            bool output = FALSE;
            if (data_array_[a] > data_array_[b]){
                output = TRUE;
            } else if (data_array_[a] == data_array_[b]){
                if (a > b){
                    output = TRUE;
                }
            }
            return output;
        });

        return buffer;
    }
};


#endif //MIGRA_GAPBS_DATA_MAP_H
