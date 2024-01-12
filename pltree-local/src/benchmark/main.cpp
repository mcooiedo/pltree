#include <iostream>

#include "pltree.h"
#include "flags.h"
#include "bench_util.h"
#include "sosd_util.h"
#include "operation.h"
#define CHECK 0

using Alloc=my_alloc::MemPoolAllocator;
Alloc::MetaType * Alloc::meta_= nullptr;
size_t Alloc::piece_size_[MemPoolAllocator::POOL_CNT]= {0};
size_t Alloc::cur_blk_[MemPoolAllocator::POOL_CNT]= {0};
size_t Alloc::max_blk_[MemPoolAllocator::POOL_CNT]= {0};
bool Alloc::recover=false;
Spinlock Alloc::alloc_mtx[MemPoolAllocator::POOL_CNT+1];
PMEMobjpool* Alloc::_pop= nullptr;
char* Alloc::buff_[MemPoolAllocator::POOL_CNT][MemPoolAllocator::PIECE_CNT] = {};
char* Alloc::buff_aligned_[MemPoolAllocator::POOL_CNT][MemPoolAllocator::PIECE_CNT] = {};
using  namespace plt;
// Modify these if running your own workload
#define PAYLOAD_TYPE char *
// Global parameters
std::string keys_file_path;
std::string keys_file_type;
std::string keys_type;
int init_num_keys;
int workload_keys; // Number of keys to operation in the workload
int total_num_keys;
std::string operation;
double insert_frac;
std::string lookup_distribution;
bool positive_search = true;
bool using_epoch = false;
bool skip_bulkload = false;
bool random_shuffle_ = false;
bool sort_bulkload = true;
int thread_num;

double theta;
int batch_size = 10000000;
template <typename T, typename P>
PLTree<T,P> *generate_index()
{
    PLTree<T,P> *index = nullptr;
    return index;
}
// Run function, to select the workload to run
// Benchmark
template <class T, class P>
void GeneralBench(PLTree<T,P> *index, Range *rarray, int thread_num, void (*test_func)(PLTree<T,P> *, struct Range *))
{
    std::thread *thread_array[1024];
    double duration;
    finished = false;
    bar_a = 1;
    bar_b = thread_num;
    bar_c = thread_num;



    for (uint64_t i = 0; i < thread_num; ++i)
    {
        thread_array[i] = new std::thread(*test_func, index, &rarray[i]);
    }

    while (LOAD(&bar_b) != 0)
        ;                                    // Spin
    std::unique_lock<std::mutex> lck(mtx); // get the lock of condition variable

    STORE(&bar_a, 0); // start test
    while (!finished)
    {
        cv.wait(lck); // go to sleep and wait for the wake-up from child threads
    }

    for (int i = 0; i < thread_num; ++i)
    {
        thread_array[i]->join();
        delete thread_array[i];
    }

    //});

    double total_throughput = 0;
    double longest_time = 0;
    uint64_t total_num = 0;
    for (int i = 0; i < thread_num; i++)
    {
        if (longest_time < rarray[i].total_time)
        {
            longest_time = rarray[i].total_time;
        }
        total_num += rarray[i].num;
    }

    total_throughput += total_num / longest_time * 1e9;

    std::cout << "\tcumulative throughput:\t"
              << total_throughput << " ops/sec"
              << std::endl;

}
double mean(std::vector<double> arr){
    double sum=0;
    for(auto item :arr){
        sum+=item;
    }
    return sum/arr.size();
}
std::pair<double ,double> getmimmax(std::vector<double> arr){
    std::sort(arr.begin(),arr.end());
    std::pair<double ,double> p;
    p.first=arr[0];
    p.second=arr.back();
    return p;
}
template <typename T,typename P=PAYLOAD_TYPE>
void Run()
{
    typedef struct KV<T,P> KV;
    // Read keys from file
    T *keys = new T[total_num_keys];
    if (keys_file_type == "binary")
    {
        load_binary_data(keys, total_num_keys, keys_file_path);
    }
    else if (keys_file_type == "text")
    {
        load_text_data(keys, total_num_keys, keys_file_path);
    }
    else if (keys_file_type == "sosd")
    {
        // Benchmark on SOSD data, using SOSD's loading method
        std::vector<T> my_keys = util::load_data<T>(keys_file_path);
        bool unique_keys_ = util::is_unique<T>(my_keys);
        if (unique_keys_)
            std::cout << "data is unique" << std::endl;
        else
            std::cout << "data contains duplicates" << std::endl;
        T *copy_keys = &my_keys[0];
        memcpy(keys, copy_keys, sizeof(T) * total_num_keys);
        random_shuffle_ = true;
    }
    else
    {
        std::cerr << "--keys_file_type must be either 'binary' or 'text'"
                  << std::endl;
        return;
    }

    if (random_shuffle_)
    {
        std::random_shuffle(&keys[0], &keys[total_num_keys - 1]);
    }

    // Combine bulk loaded keys with generated values
    std::vector<KV> kvs;
    std::mt19937_64 gen_payload(std::random_device{}());


   // PLTree<T,P> *index = generate_index<T, P>();
    PLTree<T,P> *index;
    my_alloc::MemPoolAllocator::createMempool(layout_name,pool_name,pool_size);
   // exit(-1);

    void* index_pm=Alloc::GetRoot(sizeof(PLTree<T,P>));

    if(Alloc::isRecover()&&CHECK){
        skip_bulkload=true;
    }
    new(index_pm)PLTree<T,P>(CHECK);
    Alloc::Persist(&index_pm,sizeof(PLTree<T,P>));
    // index=index_pm;
    index = reinterpret_cast<PLTree<T,P> *>(index_pm);
    char *rammem;
#ifdef LARGE_DRAM
    //my_alloc::allocMem2(rammem,allocate_mem);
    rammem = new char[allocate_mem+63];
    // void* ptr = new char[size + 63];
    auto aligned_ptr = reinterpret_cast<char*>((reinterpret_cast<uintptr_t>(rammem) + 63) & ~63);
    start_mem = (char *)aligned_ptr;
    curr_mem = start_mem;
    thread_mem_start_addr = (char *)aligned_ptr + MEM_OF_MAIN_THREAD;
#endif
    auto initKey=new T[init_num_keys];
    memcpy(initKey,keys,init_num_keys*sizeof(T));
    // Bulk loading keys to the index
    if (sort_bulkload)
    {
        multi_threaded_sort(initKey, init_num_keys,110);/** */
        std::sort(initKey,initKey+init_num_keys);
        /*std::sort(kvs.begin(), kvs.begin() + init_num_keys,
                   [](auto const &a, auto const &b)
                   { return a.first < b.first; });*/
    }
    for (int i = 0; i < init_num_keys; i++)
    {
        kvs.push_back(KV{initKey[i],reinterpret_cast<P>(&initKey[i])});
    }
delete[] initKey;
//    std::vector<T> init_keys;
//    std::vector<P> init_values;
//    for(auto kv:kvs){
//        init_keys.push_back(kv.first);
//        init_values.push_back(kv.second);
//    }

    if (!skip_bulkload)
    {
        std::cout << "Start the bulk load" << std::endl;
        std::cout << "The min key = " << kvs[0].first << std::endl;
        std::cout << "The max key = " << kvs[init_num_keys - 1].first<< std::endl;
        auto start=curr_mem;
        auto workload_start_time = std::chrono::high_resolution_clock::now();
        index->BulkLoad(kvs);
        auto workload_end_time = std::chrono::high_resolution_clock::now();
        kvs.clear();
        kvs.shrink_to_fit();
        std::cout << "End the bulk load:" <<std::chrono::duration_cast<std::chrono::nanoseconds>(workload_end_time -
                                                                                                 workload_start_time).count()<<"nanoseconds" <<std::endl;
        std::cout<<"memcost="<<(curr_mem-start)/1024.0/1024/1024<<std::endl;

    }

    int i = init_num_keys;
    int num_pre_insertion = i;

    // Generate the workload
    // Mixed workload (mixing search and insert), search/erase/update workload (existence keys), insert workload (new keys)
    uint64_t generate_num = workload_keys;
    uint64_t total_num=total_num_keys;
    char *workload = generate_workload(operation, generate_num, keys, num_pre_insertion, insert_frac, batch_size, lookup_distribution, theta,total_num,positive_search);
    int num_cores = std::thread::hardware_concurrency();
    // Partition the workload
    Range *range_array = new Range[thread_num];
    if (operation == "mixed")
    {
        // Benchmark mixed operationwq2
        typedef std::pair<int, T> OPT;
        OPT *my_workload = reinterpret_cast<OPT *>(workload);
        auto partition_num = generate_num / thread_num;
        for (int i = 0; i < (thread_num - 1); ++i)
        {
            range_array[i].id = i*2;//2 * i;
            range_array[i].num = partition_num;
            range_array[i].workload = reinterpret_cast<char *>(my_workload + i * partition_num);
        }

        range_array[thread_num - 1].id = 2 * (thread_num - 1);//
        range_array[thread_num - 1].num = generate_num - partition_num * (thread_num - 1);
        range_array[thread_num - 1].workload = reinterpret_cast<char *>(my_workload + (thread_num - 1) * partition_num);
    }
    else
    {
        // Benchmark single operation
        T *my_workload = reinterpret_cast<T *>(workload);
        auto partition_num = generate_num / thread_num;
        for (int i = 0; i < (thread_num - 1); ++i)
        {
            range_array[i].id = 2 *i;//
            range_array[i].num = partition_num;
            range_array[i].workload = reinterpret_cast<char *>(my_workload + i * partition_num);
        }

        range_array[thread_num - 1].id = 2 * (thread_num - 1);//
        range_array[thread_num - 1].num = generate_num - partition_num * (thread_num - 1);
        range_array[thread_num - 1].workload = reinterpret_cast<char *>(my_workload + (thread_num - 1) * partition_num);
    }

    // Run the benchmark
    if (using_epoch == true)
    {
       /* if (operation == "mixed")
        {
            GeneralBench<T, PAYLOAD_TYPE>(index, range_array, thread_num, &concurr_mixed);
        }
        else if (operation == "insert")
        {
            GeneralBench<T, PAYLOAD_TYPE>(index, range_array, thread_num, &concurr_insert);
        }
        else if (operation == "search")
        {
            GeneralBench<T, PAYLOAD_TYPE>(index, range_array, thread_num, &concurr_search);
        }
        else if (operation == "erase")
        {
            GeneralBench<T, PAYLOAD_TYPE>(index, range_array, thread_num, &concurr_erase);
        }
        else if (operation == "update")
        {
            GeneralBench<T, PAYLOAD_TYPE>(index, range_array, thread_num, &concurr_update);
        }
        else if (operation == "range")
        {
            GeneralBench<T, PAYLOAD_TYPE>(index, range_array, thread_num, &concurr_range);
        }
        else
        {
            std::cout << "Unknown operation " << operation << std::endl;
        }*/
    }
    else
    {
        if (operation == "mixed")
        {
            GeneralBench<T, PAYLOAD_TYPE>(index, range_array, thread_num, &mixed_without_epoch);
        }
        else if (operation == "insert")
        {
            GeneralBench<T, PAYLOAD_TYPE>(index, range_array, thread_num, &insert_without_epoch);
        }
        else if (operation == "search")
        {
            GeneralBench<T, PAYLOAD_TYPE>(index, range_array, thread_num, &search_without_epoch);
        }
        else if (operation == "erase")
        {
            GeneralBench<T, PAYLOAD_TYPE>(index, range_array, thread_num, &erase_without_epoch);
        }
        else if (operation == "update")
        {
            GeneralBench<T, PAYLOAD_TYPE>(index, range_array, thread_num, &update_without_epoch);
        }
        else if (operation == "range")
        {
            GeneralBench<T, PAYLOAD_TYPE>(index, range_array, thread_num, &range_without_epoch);
        }
        else
        {
            std::cout << "Unknown operation " << operation << std::endl;
        }
    }
    //auto p=getmimmax(my_alloc::statistic.vars);
/*    std::cout << "initAlloc = " << my_alloc::statistic.initAlloc <<std::endl;
    std::cout << "initAllocSize = " << my_alloc::statistic.initAllocSize <<std::endl;
    std::cout << "my_alloc::statistic.Inserttobuffer = " << my_alloc::statistic.Inserttobuffer <<std::endl;
    std::cout << "segNUM = " << my_alloc::statistic.segNUM <<std::endl;
    std::cout << "persisnum = " << my_alloc::statistic.persisnum <<std::endl;
    std::cout << "allocNUM = " << my_alloc::statistic.allocnum <<std::endl;
    std::cout << "maxlevel= " << my_alloc::statistic.maxlevel<<std::endl;
    std::cout << "treeHeight  = " << my_alloc::statistic.maxdegree << " ;maxoverflow = " << my_alloc::statistic.maxoverflow <<std::endl;
    std::cout << "nodeused = " << my_alloc::statistic.bufferused <<std::endl;
    std::cout << "expand = " << my_alloc::statistic.expand << ";newnode = " << my_alloc::statistic.newnode <<std::endl;
    std::cout << " blockvisit= " <<  my_alloc::statistic.bufferused<< ";buffertosgement = " << my_alloc::statistic.buffertosgement<<std::endl;
    std::cout <<  "=0= " << my_alloc::statistic.block0<< ";<15 = " << my_alloc::statistic.block15<<";>15 = " << my_alloc::statistic.blockmore15<<std::endl;*/
   // std::cout << "mean = " << mean(my_alloc::statistic.vars)<<"min = " << p.first<<"max = " << p.second<<std::endl;
    //std::cout<<"rmse = "<<sqrt(my_alloc::statistic.loss/num)<<std::endl;
    //std::cout<<"Height = "<<my_alloc::statistic.height<<std::endl;
    delete[] keys;
    index->freeSegMetas();
    index->closePool();
#ifdef LARGE_DRAM
    //delete[] rammem;
    //exit(-1);
#endif
}

int main(int argc, char *argv[]) {

    //set_affinity(0);

    // Get the flag from user
    auto flags = parse_flags(argc, argv);
    keys_file_path = get_required(flags, "keys_file");
    keys_file_type = get_required(flags, "keys_file_type");
    keys_type = get_required(flags, "keys_type");
    init_num_keys = stoi(get_required(flags, "init_num_keys"));
    workload_keys = stoi(get_required(flags, "workload_keys")); // Number of operations in the workload
    total_num_keys = stoi(get_required(flags, "total_num_keys"));
    operation = get_required(flags, "operation"); // Which operation to evalaute
    insert_frac = stod(get_with_default(flags, "insert_frac", "0.5"));
    theta = stod(get_with_default(flags, "theta", "0.99"));

    lookup_distribution =
            get_with_default(flags, "lookup_distribution", "zipf");
    int pos_search= stoi(get_required(flags, "positive_search"));
    if (pos_search)
    {
        positive_search = true;
    }
    else
    {
        positive_search = false;
    }

    random_shuffle_ = get_boolean_flag(flags, "random_shuffle");
    int sort_flag = stoi(get_required(flags, "sort_bulkload"));
    if (sort_flag)
    {
        sort_bulkload = true;
    }
    else
    {
        sort_bulkload = false;
    }
    skip_bulkload = get_boolean_flag(flags, "skip_bulkload");
    thread_num = stoi(get_with_default(flags, "thread_num", "1"));

    // Print some critical information
    std::cout << "The key type is " << keys_type << std::endl;
   if (keys_type == "double")
    {
        Run<uint64_t>();
    }
    else
    {
        Run<uint64_t>(); // Other keys are regarded as uint64_t
    }
    return 0;


}
