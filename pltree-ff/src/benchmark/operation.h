// Define the benchmarking operations
#pragma once

#include "pltree.h"
#include <cstdint>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <sched.h>
#include <chrono>
#include "../util/utils.h"
#include "../util/allocator.h"
#include <cmath>

#define EPOCH_DURATION 1000
#define TO_SCAN 100
using namespace plt;
/*** Benchmark Control Function ***/
char *thread_space_start_addr;
__thread char *start_addr;
__thread char *curr_addr;
char *thread_mem_start_addr;
__thread char *start_mem;
__thread char *curr_mem;
const uint64_t SPACE_PER_THREAD = 2 * 1024ULL * 1024ULL*1024ULL;
const uint64_t SPACE_OF_MAIN_THREAD = 10ULL * 1024ULL * 1024ULL * 1024ULL;
const uint64_t MEM_PER_THREAD = 512ULL * 1024ULL * 1024ULL;
const uint64_t MEM_OF_MAIN_THREAD = 1ULL * 1024ULL * 1024ULL * 1024ULL;

bool finished = false;
int bar_a, bar_b, bar_c;
std::mutex mtx;
std::condition_variable cv;

struct operation_record_t {
  uint64_t number = 0;
  uint64_t dummy[8]; /*patch to a cacheline size, avoid false sharing*/
};

operation_record_t operation_record[1024]; // Used for sampling

//void set_affinity(uint32_t idx) {
//  cpu_set_t my_set;
//  CPU_ZERO(&my_set);
//  CPU_SET(idx, &my_set);
//  sched_setaffinity(0, sizeof(cpu_set_t), &my_set);
//}
void set_addr(int workerid){
#ifdef LARGE_DRAM
    start_mem = thread_mem_start_addr + workerid * MEM_PER_THREAD;
    curr_mem = start_mem;
#endif
}
void set_affinity(uint32_t idx) {
   // idx=idx/2;
    cpu_set_t my_set;
    CPU_ZERO(&my_set);
    if (idx < 56) {
        CPU_SET(idx, &my_set);
    } else {
        CPU_SET(idx, &my_set);
    }
    pthread_t current_thread = pthread_self();
    pthread_setaffinity_np(current_thread, sizeof(cpu_set_t), &my_set);
    //sched_setaffinity(0, sizeof(cpu_set_t), &my_set);
}
inline void spin_wait() {
  SUB(&bar_b, 1);
  while (LOAD(&bar_a) == 1)
    ; /*spinning*/
}

inline void end_notify() {
  if (SUB(&bar_c, 1) == 0) {
    std::unique_lock<std::mutex> lck(mtx);
    finished = true;
    cv.notify_one();
  }
}

struct Range {
  int id; // To specify which core to attach
  uint64_t num; // Total insertion num
  char *workload; // Real workload, and please ensure this is the start of 
  double total_time; // Consumed time
};

/*** Benchmark Single Operation ***/

// Without Epoch operation
template <class T, class P>
void insert_without_epoch(PLTree<T,P> *index, struct Range *_range) {
    set_affinity(_range->id);
    uint64_t num = _range->num;
    T *key_array = reinterpret_cast<T*>(_range->workload);
    int fail_insert = 0;

    spin_wait();
    auto workload_start_time = std::chrono::high_resolution_clock::now();

    for (uint64_t i = 0; i < num; ++i) {
      bool ret = index->insert(key_array[i], reinterpret_cast<P>(&key_array[i]));
     // index->
      if(!ret) fail_insert++;
    }
    
    auto workload_end_time = std::chrono::high_resolution_clock::now();
    _range->total_time =
      std::chrono::duration_cast<std::chrono::nanoseconds>(workload_end_time -
                                                            workload_start_time)
          .count();
    std::cout << "Fail insert = " << fail_insert << std::endl;
    end_notify();
}


// Without Epoch operation
template <class T, class P>
void search_without_epoch(PLTree<T,P> *index, struct Range *_range) {
    set_addr(_range->id/2);
    set_affinity(_range->id);
    uint64_t num = _range->num;
    T *key_array = reinterpret_cast<T *>(_range->workload);
    uint64_t found = 0;
    P payload;
    uint64_t total_num = 0;

    spin_wait();
    auto workload_start_time = std::chrono::high_resolution_clock::now();

    for (uint64_t i = 0; i < num; ++i) {
      bool ret = index->search(key_array[i], payload);
      if(ret){
        found++;
        total_num += reinterpret_cast<uint64_t>(payload);
      }
    }
    
    auto workload_end_time = std::chrono::high_resolution_clock::now();
    _range->total_time =
      std::chrono::duration_cast<std::chrono::nanoseconds>(workload_end_time -
                                                            workload_start_time)
          .count();
    std::cout << "Not found during search = " << num - found << " with total " << total_num << std::endl;   
    end_notify();
}


// Without Epoch operation
template <class T, class P>
void range_without_epoch(PLTree<T,P> *index, struct Range *_range) {
    set_addr(_range->id/2);
    set_affinity(_range->id);
    uint64_t num = _range->num;
    T *key_array = reinterpret_cast<T *>(_range->workload);
    uint64_t found = 0;
    typedef struct KV<T,P> V;
    static thread_local std::map<T, P> results;
    uint64_t not_enough = 0;

    spin_wait();
    auto workload_start_time = std::chrono::high_resolution_clock::now();

    for (uint64_t i = 0; i < num; ++i) {
        int deleted= 0;
      found = index->rangeScan(key_array[i], TO_SCAN, results);
      if(found != TO_SCAN) {
        not_enough++;
      }
    }
    
    auto workload_end_time = std::chrono::high_resolution_clock::now();
    _range->total_time =
      std::chrono::duration_cast<std::chrono::nanoseconds>(workload_end_time -
                                                            workload_start_time)
          .count();
    std::cout << "Not enough records during search = " << not_enough << std::endl;   
    end_notify();
}




// Without Epoch operation
template <class T, class P>
void erase_without_epoch(PLTree<T,P> *index, struct Range *_range) {
    set_addr(_range->id/2);
    set_affinity(_range->id);
    uint64_t num = _range->num;
    T *key_array = reinterpret_cast<T *>(_range->workload);
    uint64_t erased = 0;

    spin_wait();
    auto workload_start_time = std::chrono::high_resolution_clock::now();

    for (uint64_t i = 0; i < num; ++i) {
      if(index->erase(key_array[i])){
        erased++;
      }
    }
    auto workload_end_time = std::chrono::high_resolution_clock::now();
    _range->total_time =
      std::chrono::duration_cast<std::chrono::nanoseconds>(workload_end_time -
                                                            workload_start_time)
          .count();
    std::cout << "Not found during erase = " << num - erased << std::endl;
    end_notify();
}


template <class T, class P>
void update_without_epoch(PLTree<T,P> *index, struct Range *_range) {
    set_addr(_range->id/2);
    set_affinity(_range->id);
    uint64_t num = _range->num;
    T *key_array = reinterpret_cast<T *>(_range->workload);
    uint64_t updates = 0;

    spin_wait();
    auto workload_start_time = std::chrono::high_resolution_clock::now();

    for (uint64_t i = 0; i < num; ++i) {
      if(index->update(key_array[i], reinterpret_cast<P>(&key_array[i]))){
        updates++;
      }
    }
    
    auto workload_end_time = std::chrono::high_resolution_clock::now();
    _range->total_time =
      std::chrono::duration_cast<std::chrono::nanoseconds>(workload_end_time -
                                                            workload_start_time)
          .count();
    std::cout << "Not found during update = " << num - updates << std::endl;
    end_notify();
}

/*** Benchmark Mixed Operation ***/

// Mixed Insert and Search operations

template <class T, class P>
void mixed_without_epoch(PLTree<T,P> *index, struct Range *_range) {
    set_addr(_range->id/2);
    set_affinity(_range->id);
    uint64_t num = _range->num;
    typedef std::pair<int, T> OPT;
    OPT *key_array = reinterpret_cast<OPT *>(_range->workload);
    uint64_t found = 0;
    thread_local P payload;
    spin_wait();
    auto start1=curr_mem;
    auto workload_start_time = std::chrono::high_resolution_clock::now();
    for (uint64_t i = 0; i < num; ++i) {
        auto op = key_array[i].first;
        auto key = key_array[i].second;
        if (op == 0) {
           bool ret = index->insert(key, reinterpret_cast<P>(&key));
        } else {
         found += index->search(key, payload);

        }
    }
    
    auto workload_end_time = std::chrono::high_resolution_clock::now();
    _range->total_time =
      std::chrono::duration_cast<std::chrono::nanoseconds>(workload_end_time -
                                                            workload_start_time).count();
  /*  std::vector<double> arr;
    auto seg=index->GetSegmentHead();
    auto tail=index->GetSegmentTail();
    int total=0;
    while (seg!=NULL){

        for(int i=0;i<seg->NumArray();++i){
            arr.push_back(seg->GetInsertNumAt(i));
            total+=seg->GetInsertNumAt(i);
        }
        seg=seg->next_segment();
    }
    double aver=(double)total/my_alloc::statistic.find;
    double t=0;
    //std::sort(arr.begin(),arr.end());
    for(int i=0;i<arr.size();++i){
        if(arr[i]!=0)
        t+=(arr[i]-aver)*(arr[i]-aver);
    }
    double avert=t/my_alloc::statistic.find;
    std::cout << "err = " << avert<< std::endl;
    std::cout << " my_alloc::statistic.find= " <<  my_alloc::statistic.find << std::endl;
*/

    std::cout<<"memcost="<<(curr_mem-start1)/1024.0/1024/1024<<std::endl;
    std::cout << "Found during mixed = " << found << std::endl;
   // std::cout<<"rmse = "<<sqrt(my_alloc::statistic.loss/num)<<std::endl;
    end_notify();
}