//
// Created by mcooiedo on 22-12-25.
//

#ifndef PLTREE_BENCH_UTIL_H
#define PLTREE_BENCH_UTIL_H

#include "zipf.h"
#include "general_zipf.h"
#include "../util/uniform.h"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <future>
template <typename IntegerType>
void multi_threaded_sort(IntegerType* in_data, size_t data_len, size_t thread_num)
{
    // compile time type check


    // under such conditions, multi-thread makes no sense
    // call std::sort directly
    if(data_len <= 1 || thread_num == 0
       || data_len < (thread_num+1)*(thread_num+1))
    {
        std::sort(in_data, in_data+data_len);
        return;
    }

    /* one thread sort one chunk
     * main thread sort the last chunk */
    size_t chunk_size = data_len/(thread_num+1);
    if(data_len%(thread_num+1) != 0)
        ++chunk_size;

    // for threads synchronize
    auto sort_promise = new std::promise<void>[thread_num];
    auto sort_future = new std::future<void>[thread_num];
    for(int i=0; i<thread_num; ++i)
        sort_future[i] = sort_promise[i].get_future();

    // create threads
    for(size_t i=0; i<thread_num; ++i){
        std::thread th([=]{
            std::sort(in_data + i*chunk_size, in_data + (i+1)*chunk_size);
            sort_promise[i].set_value();
        });
        th.detach();
    }

    // sort the last chunk
    std::sort(in_data + chunk_size*thread_num, in_data + data_len);

    // before wait and block, do things not based on data
    auto out_data = new IntegerType[data_len];
    auto index = new size_t[thread_num + 1];
    for (int i=0; i<thread_num + 1; ++i)
        index[i] = i * chunk_size;

    // wait for all threads
    for(size_t i=0; i<thread_num; ++i)
        sort_future[i].wait();

    delete[] sort_future;
    delete[] sort_promise;

    // do merge sort
    for(size_t i = 0; i < data_len; ++i)
    {
        IntegerType min_index;
        IntegerType min_num = std::numeric_limits<IntegerType>::max();

        // traverse every chunk and find the minimum
        for(size_t j=0; j<thread_num; ++j)
        {
            if((index[j] < (j+1)*chunk_size)
               && (in_data[index[j]] < min_num))
            {
                min_index = j;
                min_num = in_data[index[j]];
            }
        }
        if(index[thread_num] < data_len
           && (in_data[index[thread_num]] < min_num))
        {
            min_index = thread_num;
        }

        out_data[i] = in_data[index[min_index]];
        index[min_index]++;
    }

    std::copy(out_data, out_data + data_len, in_data);
    delete[] out_data;
}

template <class T>
bool load_binary_data(T data[], int length, const std::string& file_path) {
    std::ifstream is(file_path.c_str(), std::ios::binary | std::ios::in);
    if (!is.is_open()) {
        return false;
    }
    is.read(reinterpret_cast<char*>(data), std::streamsize(length * sizeof(T)));
    is.close();
    return true;
}

template <class T>
bool load_text_data(T array[], int length, const std::string& file_path) {
    std::ifstream is(file_path.c_str());
    if (!is.is_open()) {
        return false;
    }
    int i = 0;
    std::string str;
    while (std::getline(is, str) && i < length) {
        std::istringstream ss(str);
        ss >> array[i];
        i++;
    }
    is.close();
    return true;
}

template <class T>
T* get_search_keys(T array[], int num_keys, int num_searches,bool positive_search,int total_key_num) {
    std::mt19937_64 gen(std::random_device{}());
    std::uniform_int_distribution<int> dis(0, positive_search?(num_keys - 1):(total_key_num-1));
    auto* keys = new T[num_searches];
    for (int i = 0; i < num_searches; i++) {
        int pos = dis(gen);
        keys[i] = array[pos];
    }
    return keys;
}
template <class T>
T* get_search_keys(T array[], int num_keys, int num_searches) {
    std::mt19937_64 gen(std::random_device{}());
    std::uniform_int_distribution<int> dis(0, num_keys - 1);
    auto* keys = new T[num_searches];
    for (int i = 0; i < num_searches; i++) {
        int pos = dis(gen);
        keys[i] = array[pos];
    }
    return keys;
}
template <class T>
T* get_search_keys_zipf(T array[], int num_keys, int num_searches) {
    auto* keys = new T[num_searches];
    ScrambledZipfianGenerator zipf_gen(num_keys);
    for (int i = 0; i < num_searches; i++) {
        int pos = zipf_gen.nextValue();
        keys[i] = array[pos];
    }
    return keys;
}

template <class T>
T* get_search_keys_zipf_with_theta(T array[], int num_keys, int num_searches, double theta) {
    auto* keys = new T[num_searches];
    std::default_random_engine generator_;
    zipfian_int_distribution<int> dis(0, num_keys - 1, theta);
    for (int i = 0; i < num_searches; i++) {
        int pos = dis(generator_);
        keys[i] = array[pos];
    }
    std::random_shuffle(&keys[0], &keys[num_keys - 1]);
    return keys;
}

template <class T>
char* generate_workload(std::string& operation_, uint64_t& generate_num, T* keys, uint64_t base_insert_num, double insert_frac_,
                        int num_lookups_per_batch, const std::string& distribution = "uniform", double theta = 0.99, uint64_t total_key_num=200000000, bool postive_search=true){
    char *workload = nullptr;

    if(operation_ == "mixed"){
        typedef std::pair<int, T> OP;
        unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
        std::mt19937_64 g2 (seed1);
        uint64_t u64Random = g2();
        UniformRandom rng(u64Random);
        uint32_t random;
        uint32_t insert_sign = static_cast<uint32_t>(insert_frac_ * 100);

        uint64_t my_inserts = 0;
        uint64_t my_searches = 0;
        uint64_t mixed_generate_num = 0;
        mixed_generate_num = generate_num;
        // generate_num only means the approximate number of inserts in mixed workload
        if(insert_sign != 0){
            mixed_generate_num = static_cast<uint64_t>(generate_num / insert_frac_);
        }
        auto mixed_workload = new OP[mixed_generate_num];
        T* lookup_keys = nullptr;
        if(num_lookups_per_batch != 0){
            lookup_keys = get_search_keys(keys, base_insert_num, num_lookups_per_batch,postive_search,total_key_num);
        }

        for(uint64_t i = 0; i < mixed_generate_num; ++i){
            random = rng.next_uint32() % 100;
            if (random < insert_sign) { /*insert*/
                mixed_workload[i].first = 0;
                mixed_workload[i].second = keys[(base_insert_num + my_inserts)%total_key_num];//__bug
                my_inserts++;
                if(my_inserts == generate_num) break;
            } else { /*get*/
                mixed_workload[i].first = 1;
                mixed_workload[i].second = lookup_keys[my_searches % num_lookups_per_batch];
                my_searches++;
                if (((my_searches % num_lookups_per_batch) == 0) && (num_lookups_per_batch != 0))
                {
                    delete [] lookup_keys;
                    lookup_keys = get_search_keys(keys, base_insert_num, num_lookups_per_batch,postive_search,total_key_num);
                }
            }
        }

        if(lookup_keys != nullptr){
            delete [] lookup_keys;
        }

        generate_num = my_inserts + my_searches;
        std::cout << "Generated mixed inserts = " << generate_num << ": #inserts = " << my_inserts << "; #searches =" << my_searches << std::endl;
        workload = reinterpret_cast<char*>(mixed_workload);
    } else if ((operation_ == "insert") || (operation_ == "debug") || (operation_ == "var_insert")) {
        // Reuse the existing key array
        auto my_workload = new T[generate_num];
        memcpy(my_workload, keys + base_insert_num, sizeof(T) * generate_num);
        workload = reinterpret_cast<char*>(my_workload);
        if(distribution == "zipf"){
            std::cout << "Only support zipf for update/search/range-query" << std::endl;
            exit(-1);
        }
    } else if (operation_ == "erase"){
        // Need to avoid the duplicate, and ensure these keys already existis
        auto my_workload = new T[generate_num];
        double erase_frac = generate_num / static_cast<double>(base_insert_num);
        uint32_t erase_sign = static_cast<uint32_t>(erase_frac * 10);
        uint64_t erase_idx = 0;

        for(uint64_t i = 0; i < base_insert_num; i += 10){
            if((base_insert_num - i - 10) < (generate_num - erase_idx)){
                uint64_t end = i + generate_num - erase_idx;
                for (uint64_t j = i; j < end; ++j)
                {
                    my_workload[erase_idx++] = keys[j];
                }
                break;
            }

            // For every 10 keys, delete erase_sign keys
            for (uint64_t j = 0; j < erase_sign; ++j)
            {
                my_workload[erase_idx++] = keys[i + j];
            }
        }

        if (erase_idx != generate_num)
        {
            LOG_FATAL("Generate wrong number of delete keys");
        }

        std::random_shuffle(&my_workload[0], &my_workload[generate_num - 1]);
        workload = reinterpret_cast<char*>(my_workload);
    } else if (operation_ == "range_debug"){
        auto my_workload = new T[generate_num];
        memcpy(my_workload, keys, sizeof(T) * generate_num);
        std::sort(my_workload, my_workload + generate_num,
                  [](auto const& a, auto const& b) { return a < b; });
        workload = reinterpret_cast<char*>(my_workload);
    } else {
        // Update/search/range-query operation
        T *lookup_keys = nullptr;
        if(distribution == "zipf"){
            std::cout << "generate zipf workload with theta = " << theta << std::endl;
            lookup_keys = get_search_keys_zipf_with_theta(keys, base_insert_num, generate_num, theta);
        }else{
            lookup_keys = get_search_keys(keys, base_insert_num, generate_num);
        }

        workload = reinterpret_cast<char*>(lookup_keys);
    }

    return workload;
}



#endif //PLTREE_BENCH_UTIL_H
