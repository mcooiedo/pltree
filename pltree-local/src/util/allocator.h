//
// Created by mcooiedo on 22-12-1.
//

#ifndef PLTREE_ALLOCATOR_H
#define PLTREE_ALLOCATOR_H

#include <iostream>
#include <cstddef>
#include <climits>
#include <cstdlib>
#include <new>
#include <vector>
#include <string>
#include <sys/stat.h>
#include <cassert>

#include "libpmem.h"
#include "libpmemobj.h"
//#include <garbage_list.h>
#include "spinlock.h"
#include "utils.h"

#define FLUSH
#define SERVER
#define LARGE_DRAM
// In this class, I will write a custom template allocator
// Specifically, it allocates persistent memory using PMDK interface
// Moreover, need to use static member to make all allocatoion in a single memory pool
POBJ_LAYOUT_BEGIN(pmallocator);
POBJ_LAYOUT_TOID(pmallocator, char)

POBJ_LAYOUT_END(pmallocator)
static const char *layout_name = "template_pool";
static const uint64_t pool_addr = 0x5f0000000000;
#ifdef SERVER
static const char *pool_name = "/mnt/pmem0.2/template.data";
static const uint64_t pool_size = 96UL * 1024 * 1024 * 1024UL;// 20UL * 1024 * 1024 * 1024;
static const uint64_t allocate_mem = 30ULL * 1024ULL * 1024ULL * 1024ULL;
#else
static const char *pool_name = "/home/mcooiedo/文档/code/template";///mnt/pmem0/baotong/template.data";
static const uint64_t pool_size =20UL * 1024 * 1024 * 1024UL ;// 20UL * 1024 * 1024 * 1024
static const uint64_t allocate_mem = 10ULL * 1024ULL * 1024ULL * 1024ULL;
#endif
extern __thread char *start_mem;
extern __thread char *curr_mem;

namespace my_alloc {
    static void align_alloc(void **ptr, size_t size) {
        posix_memalign(ptr, 64, size);
    }

    static void align_zalloc(void **ptr, size_t size) {
        int ret = posix_memalign(ptr, 64, size);
        if (ret != 0) {
            // 内存分配失败
            LOG_FATAL("failed to allocate memory1");
        }
        memset(*ptr, 0, size);
    }
    template<typename T>
    static void allocMem1(T *&ret, int size) {
        uint64_t allocblock=size * sizeof(T)%CACHELINE_SIZE==0?size * sizeof(T):size * sizeof(T)+(CACHELINE_SIZE-size * sizeof(T)%CACHELINE_SIZE);
        // align_zalloc(&ptr, size * sizeof(T));
        ret = reinterpret_cast<T *>(curr_mem);
        curr_mem+=allocblock;
    }
    template<typename T>
    static void allocMem2(T *&ret, int size) {
        void *ptr;
        align_zalloc(&ptr, size * sizeof(T));
        ret = reinterpret_cast<T *>(ptr);
    }

    template<typename T>
    static void allocMem(T *&ret, int size) {
#ifdef LARGE_DRAM
        allocMem1(ret, size);
#else
        allocMem2(ret, size);
#endif
       // void *ptr;
        //align_zalloc(&ptr, size * sizeof(T));
       // ret = reinterpret_cast<T *>(ptr);
    }


    struct numData {
        int initAlloc = 0;
        uint64_t initAllocSize = 0;
        int allocnum = 0;
        int persisnum = 0;
        int maxoverflow = 0;
        int retry = 0;
        int bufferused = 0;
        size_t maxlevel = 0;
        int maxdegree = 0;
        int Inserttobuffer = 0;
        int expand = 0;
        int height = 0;
        std::vector<double> vars;
        std::vector<double> vars2;
        std::vector<int> vars3;
        int newnode = 0;
        int blockmore15 = 0;
        int block15 = 0;
        int block0 = 0;
        int buffertosgement = 0;
        int segmentresize = 0;
        int insertToStash = 0;
        int segNUM = 0;
        double loss;

        //计算vector的mean
        double mean() {
            double sum = 0;
            for (int i = 0; i < vars.size(); i++) {
                sum += vars[i];
            }
            return sum / vars.size();
        }
    };

    static numData statistic;

    typedef void (*DestroyCallback)(void *callback_context, void *object);

    template<class T1, class T2>
    inline void _construct(T1 *p, const T2 &value) { new(p) T1(value); }

    template<class T>
    inline void _destroy(T *ptr) { ptr->~T(); }

    // Implement a base class that has the memory pool
    class PMPool {
    public:
        static PMEMobjpool *_pop;
        bool recover = false;

        PMPool() {
        }

        ~PMPool() {
            //ClosePool();
        }

        static void ClosePool() {
            if (_pop != nullptr) {
                pmemobj_close(_pop);
            }
        }

    private:

    };

    class PMAllocator {
    public:
        static size_t allocNum;
        static size_t dallocNum;
        static PMEMobjpool *_pop;
        static bool recover;

        static bool PMPoolInit(const char *layout_name, const char *file_name, size_t pool_size) {
            recover = false;
            if (!FileExists(file_name)) {
                pool_size = pool_size + ((pool_size & ((1 << 23) - 1)) > 0 ? (1 << 23) : 0);
                LOG("creating a new pool");
                _pop = pmemobj_create(file_name, layout_name, pool_size, CREATE_MODE_RW);
                if (!_pop) {
                    LOG_FATAL("failed to create a pool;");
                }
                std::cout << "pool opened at: " << std::hex << _pop << std::dec << std::endl;
            } else {
                LOG("opening an existing pool, and trying to map to same address");
                /* Need to open an existing persistent pool */
                recover = true;
                _pop = pmemobj_open(file_name, layout_name);
                //pm_pool_ = pmemobj_open(pool_name, layout_name);
                if (_pop == nullptr) {
                    LOG_FATAL("failed to open the pool");
                }
                std::cout << "pool opened at: " << std::hex << _pop
                          << std::dec << std::endl;
            }
            return recover;
        }

        PMAllocator() {

        }

        bool isRecover() {
            return recover;
        }

        void *GetRoot(size_t size) {
            return pmemobj_direct(pmemobj_root(_pop, size));
        }

        void static Persist(void *p, size_t size) {
            //statistic.persisnum++;
            volatile char *ptr = (char *) ((unsigned long long) p & ~(CACHELINE_SIZE - 1));
            // if(fence) mfence();
            for (; ptr < (char *) p + size; ptr += CACHELINE_SIZE) {
                flush((void *) ptr);
            }
            sfence();
            //pmemobj_persist(pool._pop, p, size);
        }

        void static Flush(void *p, size_t size) {
            //statistic.persisnum++;
            volatile char *ptr = (char *) ((unsigned long long) p & ~(CACHELINE_SIZE - 1));
            for (; ptr < (char *) p + size; ptr += CACHELINE_SIZE) {
                flush((void *) ptr);
            }
            //pmemobj_flush(pool._pop,p,size);
        }

        void Free(void *p) {
            auto ptr = pmemobj_oid(p);
            pmemobj_free(&ptr);
        }



        void PMZAllocAt(PMEMoid *tmp, size_t msize) {
            //statistic.allocnum++;
            auto ret = pmemobj_zalloc(_pop, tmp, msize, TOID_TYPE_NUM(char));
            if (ret) {
                std::cout << "Fail logging: " << ret << "; Size = " << msize << std::endl;
                LOG_FATAL("Allocate: Allocation Error in PMEMoid 1");
            }

        }

        void *PMZAlloc(size_t msize) {
            //statistic.allocnum++;
            PMEMoid tmp;
            auto ret = pmemobj_zalloc(_pop, &tmp, msize, TOID_TYPE_NUM(char));
            if (ret) {
                std::cout << "Fail logging: " << ret << "; Size = " << msize << std::endl;
                LOG_FATAL("Allocate: Allocation Error in PMEMoid 1");
            }
            void *mem = pmemobj_direct(tmp);
            assert(mem != nullptr);
            return mem;

        }

        void *GetRootPtr(size_t size) {
            return pmemobj_direct(pmemobj_root(_pop, size));
        }

        inline void *pmAlloc(size_t msize) {
            //statistic.allocnum++;
            PMEMoid tmp;
            auto ret = pmemobj_alloc(_pop, &tmp, msize, TOID_TYPE_NUM(char), NULL, NULL);
            if (ret) {
                LOG_FATAL("pmAlloc failed \n");
            }
            void *mem = pmemobj_direct(tmp);
            assert(mem != nullptr);
            return mem;
        }
        // inline
        // private:




    };

    class MemPoolAllocator {

    public:

        //static const int64_t DEFAULT_POOL_SIZE = POOLSIZE * (1024UL * 1024 * 1024); // set the default pool size to be 10GB
        static const int POOL_CNT = 100;
        static const int PIECE_CNT = 64;
        static const size_t ALIGN_SIZE = 256;
        static char *buff_[POOL_CNT][PIECE_CNT];
        static char *buff_aligned_[POOL_CNT][PIECE_CNT];
        struct MetaData{
            char *buffer[PIECE_CNT];
            size_t blk_per_piece;
            size_t cur_blk;
        };
        struct MetaType {
            MetaData metas[POOL_CNT];
            // entrance of DS in buffer
            void *entrance;
        };
        static MetaType *meta_;
        static size_t piece_size_[POOL_CNT];
        static size_t cur_blk_[POOL_CNT];
        static size_t max_blk_[POOL_CNT];
        static Spinlock alloc_mtx[POOL_CNT+1];

        bool file_exist(const char *pool_path) {
            struct stat buffer;
            return (stat(pool_path, &buffer) == 0);
        }

    public:
        static bool recover;
        static PMEMobjpool *_pop;

        /*
         *  Construct a PM allocator, map a pool file into virtual memory
         *  @param filename     pool file name
         *  @param recover      if doing recover, false for the first time
         *  @param layout_name  ID of a group of allocations (in characters), each ID corresponding to a root entry
         *  @param pool_size    pool size of the pool file, vaild if the file doesn't exist
         */

        static bool isRecover() {
            return recover;
        }

        static bool PMPoolInit(const char *layout_name, const char *file_name, size_t pool_size) {
            recover = false;
            if (!FileExists(file_name)) {
                pool_size = pool_size + ((pool_size & ((1 << 23) - 1)) > 0 ? (1 << 23) : 0);
                LOG("creating a new pool");
                _pop = pmemobj_create(file_name, layout_name, pool_size, CREATE_MODE_RW);
                if (!_pop) {
                    LOG_FATAL("failed to create a pool;");
                }
                std::cout << "pool opened at: " << std::hex << _pop << std::dec << std::endl;
            } else {
                LOG("opening an existing pool, and trying to map to same address");
                /* Need to open an existing persistent pool */
                recover = true;
                _pop = pmemobj_open(file_name, layout_name);
                //pm_pool_ = pmemobj_open(pool_name, layout_name);
                if (_pop == nullptr) {
                    LOG_FATAL("failed to open the pool");
                }
                std::cout << "pool opened at: " << std::hex << _pop
                          << std::dec << std::endl;
            }
            return recover;
        }

        static void createMempool(const char *layout_name, const char *file_name, int64_t pool_size) {
            pool_size = pool_size + ((pool_size & ((1 << 23) - 1)) > 0 ? (1 << 23) : 0);
            recover = PMPoolInit(layout_name, file_name, pool_size);
            // align to 8MB
            if (recover == false) {

                meta_ = (MetaType *) pmemobj_direct(pmemobj_root(_pop, sizeof(MetaType)));

                // maintain volatile domain
                uint64_t alloc_size = (pool_size >> 3)*6/POOL_CNT; // x of the pool is used as block alloction
                for(int j= 0; j < POOL_CNT; ++j) {
                    for (int i = 0; i < PIECE_CNT; i++) {
                        buff_[j][i] = (char *) mem_alloc(alloc_size / PIECE_CNT);
                        buff_aligned_[j][i] = (char *) ((uint64_t) buff_[j][i] +
                                                     ((uint64_t) buff_[j][i] % ALIGN_SIZE == 0 ? 0 : (ALIGN_SIZE -
                                                                                                   (uint64_t) buff_[j][i] %
                                                                                                   ALIGN_SIZE)));
                    }
                    piece_size_ [j]= (alloc_size / PIECE_CNT) / ALIGN_SIZE - 1;
                    max_blk_[j] = piece_size_ [j]* PIECE_CNT;
                    cur_blk_[j] = 0;
                    // initialize meta_
                    for (int i = 0; i < PIECE_CNT; i++)
                        meta_->metas[j].buffer[i] = relative(buff_[j][i]);
                    meta_->metas[j].blk_per_piece = piece_size_[j];
                    meta_->metas[j].cur_blk = 0;
                }

                meta_->entrance = NULL;
                Flush(meta_, sizeof(MetaType));
            } else {
                meta_ = (MetaType *) pmemobj_direct(pmemobj_root(_pop, sizeof(MetaType)));
                for(int j= 0; j < POOL_CNT; ++j) {
                    for (int i = 0; i < PIECE_CNT; i++) {
                        buff_[j][i] = absolute(meta_->metas[j].buffer[i]);
                        buff_aligned_[j][i] = (char *) ((uint64_t) buff_[j][i] +
                                                     ((uint64_t) buff_[j][i] % ALIGN_SIZE == 0 ? 0 : (ALIGN_SIZE -
                                                                                                   (uint64_t) buff_[j][i] %
                                                                                                   ALIGN_SIZE)));
                    }
                    cur_blk_[j] = meta_->metas[j].cur_blk;
                    piece_size_[j] = meta_->metas[j].blk_per_piece;
                    max_blk_[j] = piece_size_[j] * PIECE_CNT;
                }

            }
        }

        MemPoolAllocator() {

        }

        ~MemPoolAllocator() {
            // pmemobj_close(pool._pop);
        }

    public:
        static void Persist(void *p, size_t size) {
            //statistic.persisnum++;
#ifdef FLUSH
            volatile char *ptr = (char *) ((unsigned long long) p & ~(CACHELINE_SIZE - 1));
            // if(fence) mfence();
            for (; ptr < (char *) p + size; ptr += CACHELINE_SIZE) {
                flush((void *) ptr);
            }
#endif
            sfence();
            //pmemobj_persist(pool._pop, p, size);
        }

        static void Flush(void *p, size_t size) {
            //statistic.persisnum++;
#ifdef FLUSH
            volatile char *ptr = (char *) ((unsigned long long) p & ~(CACHELINE_SIZE - 1));
            for (; ptr < (char *) p + size; ptr += CACHELINE_SIZE) {
                flush((void *) ptr);
            }
#endif
            //pmemobj_flush(pool._pop,p,size);
        }

        static inline void do_flush_with_double_fence(void *data, int len, bool fence = true) {
            //statistic.persisnum++;
#ifdef FLUSH
            volatile char *ptr = (char *) ((unsigned long long) data & ~(CACHELINE_SIZE - 1));
            // if(fence) mfence();
            for (; ptr < (char *) data + len; ptr += CACHELINE_SIZE) {
                flush((void *) ptr);
            }
            if (fence) sfence();
#endif           //DOFLUSH
        }

        /*
         *  Get/allocate the root entry of the allocator.
         *
         *  The root entry is the entrance of one group of allocation, each group is
         *  identified by the layout_name when constructing it.
         *
         *  Each group of allocations is a independent, self-contained in-memory structure in the pool
         *  such as b-tree or link-list
         */

        static void *GetRoot(size_t nsize) { // the root of DS stored in buff_ is recorded at meta_->entrance
            if (meta_->entrance == NULL) {
                meta_->entrance = pmAlloc(nsize);
                Flush(&meta_->entrance, sizeof(void *));
            }
            return meta_->entrance;
        }

        void PMZAllocAt(PMEMoid *tmp, size_t msize) {
            //statistic.allocnum++;
            auto ret = pmemobj_zalloc(_pop, tmp, msize, TOID_TYPE_NUM(char));
            if (ret) {
                std::cout << "Fail logging: " << ret << "; Size = " << msize << std::endl;
                LOG_FATAL("Allocate: Allocation Error in PMEMoid 1");
            }

        }

        static void *PMZAlloc(size_t nsize) {
            return pmAlloc(nsize);
        }

        /*
         *  Allocate a non-root piece of persistent memory from the mapped pool
         *  return the virtual memory address
         */
        static void *pmAlloc(size_t nsize) {

            if (nsize >= (1 << 23)) { // large than 8KB
                void *mem = mem_alloc(nsize + ALIGN_SIZE); // not aligned
                //  |  UNUSED    |HEADER|       memory you can use     |
                // mem             (mem + off)
                uint64_t offset = ALIGN_SIZE - (uint64_t) mem % ALIGN_SIZE;
                // store a header in the front
                uint64_t *header = (uint64_t *) ((uint64_t) mem + offset - 8);
                *header = offset;

                return (void *) ((uint64_t) mem + offset);
            }


            int blk_demand = (nsize + ALIGN_SIZE - 1) / ALIGN_SIZE;
            std::thread::id threadId = std::this_thread::get_id();

            auto  threadIdValue =*(unsigned int *) &threadId;
            auto poolidx = threadIdValue % POOL_CNT;
            auto startidx=poolidx;
            int full_num = 0;
            RETRY:

            if(full_num>POOL_CNT){
                printf("run out of memory\n");
                exit(-1);
            }


            if(!alloc_mtx[poolidx].trylock()){
                if (blk_demand + cur_blk_[poolidx] > max_blk_[poolidx]) {
                    if(poolidx==startidx&&full_num<POOL_CNT){
                        full_num=0;
                    }
                    full_num++;
                }
                poolidx=(poolidx+1)%POOL_CNT;
                goto RETRY;
            }
            // case 1: not enough in the buffer
            if (blk_demand + cur_blk_[poolidx] > max_blk_[poolidx]) {
                if(poolidx==startidx&&full_num<POOL_CNT){
                    full_num=0;
                }
                full_num++;
                poolidx=(poolidx+1)%POOL_CNT;
                goto RETRY;
            }
            // case 2: current piece can not accommdate this allocation
            int piece_id = cur_blk_[poolidx] / piece_size_[poolidx];
            if ((cur_blk_[poolidx] % piece_size_[poolidx] + blk_demand) > piece_size_[poolidx]) {
                void *mem = buff_aligned_[poolidx][piece_id + 1]; // allocate from a new peice
                cur_blk_[poolidx] = piece_size_[poolidx] * (piece_id + 1) + blk_demand;
                meta_->metas[poolidx].cur_blk = cur_blk_[poolidx];
                Flush(&(meta_->metas[poolidx].cur_blk), 8);

                alloc_mtx[poolidx].unlock();
                return mem;
            }
                // case 3: current piece has enough space
            else {
                void *mem = buff_aligned_[poolidx][piece_id] + ALIGN_SIZE * (cur_blk_[poolidx] % piece_size_[poolidx]);

                cur_blk_[poolidx] = cur_blk_[poolidx] + blk_demand;
                meta_->metas[poolidx].cur_blk = cur_blk_[poolidx];
                Flush(&(meta_->metas[poolidx].cur_blk), 8);
                alloc_mtx[poolidx].unlock();
                return mem;
            }
        }

        static void free(void *addr) {
/*            for (int i = 0; i < PIECE_CNT; i++) {
                uint64_t offset = (uint64_t) addr - (uint64_t) buff_aligned_[i];
                if (offset > 0 && offset < piece_size_ * ALIGN_SIZE) {
                    // the addr is in this piece, do not reclaim it
                    return;
                }
            }*/

            // larger than 4KB, reclaim it
            uint64_t *header = (uint64_t *) ((uint64_t) addr - 8);
            uint64_t offset = *header;

            alloc_mtx[POOL_CNT].lock();
            auto oid_ptr = pmemobj_oid((void *) ((uint64_t) addr - offset));
            alloc_mtx[POOL_CNT].unlock();

            TOID(char) ptr_cpy;
            TOID_ASSIGN(ptr_cpy, oid_ptr);
            POBJ_FREE(&ptr_cpy);

        }

        /*
         *  Distinguish from virtual memory address and offset in the pool
         *  Each memory piece allocated from the pool has an in-pool offset, which remains unchanged
         *  until reclaiment. We cannot ensure that the pool file is mapped at the same position at
         *  any time, so it may locate at different virtual memory addresses next time.
         *
         *  So the rule is that, using virtual memory when doing normal operations like to DRAM
         *  space, using offset to store link relationship, for exmaple, next pointer in linklist
         */

        /*
         *  convert an offset to a virtual memory address
         */
        template<typename T>
        static inline T *absolute(T *pmem_offset) {
            if (pmem_offset == NULL)
                return NULL;
            return reinterpret_cast<T *>(reinterpret_cast<uint64_t>(pmem_offset) + reinterpret_cast<char *>(_pop));
        }

        /*
         *  convert a virtual memory address to an offset
         */
        template<typename T>
        static inline T *relative(T *pmem_direct) {
            if (pmem_direct == NULL)
                return NULL;
            return reinterpret_cast<T *>(reinterpret_cast<char *>(pmem_direct) - reinterpret_cast<char *>(_pop));
        }

    private:
        static void *mem_alloc(size_t nsize) {
            //statistic.allocnum++;
            PMEMoid tmp;
            alloc_mtx[POOL_CNT].lock();
            auto ret = pmemobj_zalloc(_pop, &tmp, nsize, TOID_TYPE_NUM(char));
            if (ret) {
                std::cout << "Fail logging: " << ret << "; Size = " << nsize << std::endl;
                LOG_FATAL("Allocate: Allocation Error in PMEMoid 11");
            }
            void *mem = pmemobj_direct(tmp);
            alloc_mtx[POOL_CNT].unlock();
            assert(mem != nullptr);
            return mem;
        }
    };

    class PoolAllocator {

    public:

        //static const int64_t DEFAULT_POOL_SIZE = POOLSIZE * (1024UL * 1024 * 1024); // set the default pool size to be 10GB
        static const int POOL_CNT = 100;
        static const int PIECE_CNT = 64;
        static const size_t ALIGN_SIZE = 256;
        static char *buff_[POOL_CNT][PIECE_CNT];
        static char *buff_aligned_[POOL_CNT][PIECE_CNT];
        struct MetaData{
            char *buffer[PIECE_CNT];
            size_t blk_per_piece;
            size_t cur_blk;
        };
        struct MetaType {
            MetaData metas[POOL_CNT];
            // entrance of DS in buffer
            void *entrance;
        };
        static MetaType *meta_;
        static size_t piece_size_[POOL_CNT];
        static size_t cur_blk_[POOL_CNT];
        static size_t max_blk_[POOL_CNT];
        static Spinlock alloc_mtx[POOL_CNT+1];

        bool file_exist(const char *pool_path) {
            struct stat buffer;
            return (stat(pool_path, &buffer) == 0);
        }

    public:
        static bool recover;
        static PMEMobjpool *_pop;

        /*
         *  Construct a PM allocator, map a pool file into virtual memory
         *  @param filename     pool file name
         *  @param recover      if doing recover, false for the first time
         *  @param layout_name  ID of a group of allocations (in characters), each ID corresponding to a root entry
         *  @param pool_size    pool size of the pool file, vaild if the file doesn't exist
         */

        static bool isRecover() {
            return recover;
        }

        static bool PMPoolInit(const char *layout_name, const char *file_name, size_t pool_size) {
            recover = false;
            if (!FileExists(file_name)) {
                pool_size = pool_size + ((pool_size & ((1 << 23) - 1)) > 0 ? (1 << 23) : 0);
                LOG("creating a new pool");
                _pop = pmemobj_create(file_name, layout_name, pool_size, CREATE_MODE_RW);
                if (!_pop) {
                    LOG_FATAL("failed to create a pool;");
                }
                std::cout << "pool opened at: " << std::hex << _pop << std::dec << std::endl;
            } else {
                LOG("opening an existing pool, and trying to map to same address");
                /* Need to open an existing persistent pool */
                recover = true;
                _pop = pmemobj_open(file_name, layout_name);
                //pm_pool_ = pmemobj_open(pool_name, layout_name);
                if (_pop == nullptr) {
                    LOG_FATAL("failed to open the pool");
                }
                std::cout << "pool opened at: " << std::hex << _pop
                          << std::dec << std::endl;
            }
            return recover;
        }

        static void createMempool(const char *layout_name, const char *file_name, int64_t pool_size) {
            pool_size = pool_size + ((pool_size & ((1 << 23) - 1)) > 0 ? (1 << 23) : 0);
            recover = PMPoolInit(layout_name, file_name, pool_size);
            // align to 8MB
            if (recover == false) {

                meta_ = (MetaType *) pmemobj_direct(pmemobj_root(_pop, sizeof(MetaType)));

                // maintain volatile domain
                uint64_t alloc_size = (pool_size >> 3)*6/POOL_CNT; // x of the pool is used as block alloction
                for(int j= 0; j < POOL_CNT; ++j) {
                    for (int i = 0; i < PIECE_CNT; i++) {
                        buff_[j][i] = (char *) mem_alloc(alloc_size / PIECE_CNT);
                        buff_aligned_[j][i] = (char *) ((uint64_t) buff_[j][i] +
                                                        ((uint64_t) buff_[j][i] % ALIGN_SIZE == 0 ? 0 : (ALIGN_SIZE -
                                                                                                         (uint64_t) buff_[j][i] %
                                                                                                         ALIGN_SIZE)));
                    }
                    piece_size_ [j]= (alloc_size / PIECE_CNT) / ALIGN_SIZE - 1;
                    max_blk_[j] = piece_size_ [j]* PIECE_CNT;
                    cur_blk_[j] = 0;
                    // initialize meta_
                    for (int i = 0; i < PIECE_CNT; i++)
                        meta_->metas[j].buffer[i] = relative(buff_[j][i]);
                    meta_->metas[j].blk_per_piece = piece_size_[j];
                    meta_->metas[j].cur_blk = 0;
                }

                meta_->entrance = NULL;
                Flush(meta_, sizeof(MetaType));
            } else {
                meta_ = (MetaType *) pmemobj_direct(pmemobj_root(_pop, sizeof(MetaType)));
                for(int j= 0; j < POOL_CNT; ++j) {
                    for (int i = 0; i < PIECE_CNT; i++) {
                        buff_[j][i] = absolute(meta_->metas[j].buffer[i]);
                        buff_aligned_[j][i] = (char *) ((uint64_t) buff_[j][i] +
                                                        ((uint64_t) buff_[j][i] % ALIGN_SIZE == 0 ? 0 : (ALIGN_SIZE -
                                                                                                         (uint64_t) buff_[j][i] %
                                                                                                         ALIGN_SIZE)));
                    }
                    cur_blk_[j] = meta_->metas[j].cur_blk;
                    piece_size_[j] = meta_->metas[j].blk_per_piece;
                    max_blk_[j] = piece_size_[j] * PIECE_CNT;
                }

            }
        }

        PoolAllocator() {

        }

        ~PoolAllocator() {
            // pmemobj_close(pool._pop);
        }

    public:
        static void Persist(void *p, size_t size) {
            //statistic.persisnum++;

            volatile char *ptr = (char *) ((unsigned long long) p & ~(CACHELINE_SIZE - 1));
            // if(fence) mfence();
            for (; ptr < (char *) p + size; ptr += CACHELINE_SIZE) {
                flush((void *) ptr);
            }

            sfence();
            //pmemobj_persist(pool._pop, p, size);
        }

        static void Flush(void *p, size_t size) {
            //statistic.persisnum++;
#ifdef FLUSH
            volatile char *ptr = (char *) ((unsigned long long) p & ~(CACHELINE_SIZE - 1));
            for (; ptr < (char *) p + size; ptr += CACHELINE_SIZE) {
                flush((void *) ptr);
            }
#endif
            //pmemobj_flush(pool._pop,p,size);
        }

        static inline void do_flush_with_double_fence(void *data, int len, bool fence = true) {
            //statistic.persisnum++;
#ifdef FLUSH
            volatile char *ptr = (char *) ((unsigned long long) data & ~(CACHELINE_SIZE - 1));
            // if(fence) mfence();
            for (; ptr < (char *) data + len; ptr += CACHELINE_SIZE) {
                flush((void *) ptr);
            }
            if (fence) sfence();
#endif           //DOFLUSH
        }

        /*
         *  Get/allocate the root entry of the allocator.
         *
         *  The root entry is the entrance of one group of allocation, each group is
         *  identified by the layout_name when constructing it.
         *
         *  Each group of allocations is a independent, self-contained in-memory structure in the pool
         *  such as b-tree or link-list
         */

        static void *GetRoot(size_t nsize) { // the root of DS stored in buff_ is recorded at meta_->entrance
            if (meta_->entrance == NULL) {
                meta_->entrance = pmAlloc(nsize);
                Flush(&meta_->entrance, sizeof(void *));
            }
            return meta_->entrance;
        }

        void PMZAllocAt(PMEMoid *tmp, size_t msize) {
            //statistic.allocnum++;
            auto ret = pmemobj_zalloc(_pop, tmp, msize, TOID_TYPE_NUM(char));
            if (ret) {
                std::cout << "Fail logging: " << ret << "; Size = " << msize << std::endl;
                LOG_FATAL("Allocate: Allocation Error in PMEMoid 1");
            }

        }

        static void *PMZAlloc(size_t nsize) {
            return pmAlloc(nsize);
        }

        /*
         *  Allocate a non-root piece of persistent memory from the mapped pool
         *  return the virtual memory address
         */
        static void *pmAlloc(size_t nsize) {

            if (nsize >= (1 << 23)) { // large than 8KB
                void *mem = mem_alloc(nsize + ALIGN_SIZE); // not aligned
                //  |  UNUSED    |HEADER|       memory you can use     |
                // mem             (mem + off)
                uint64_t offset = ALIGN_SIZE - (uint64_t) mem % ALIGN_SIZE;
                // store a header in the front
                uint64_t *header = (uint64_t *) ((uint64_t) mem + offset - 8);
                *header = offset;

                return (void *) ((uint64_t) mem + offset);
            }


            int blk_demand = (nsize + ALIGN_SIZE - 1) / ALIGN_SIZE;
            std::thread::id threadId = std::this_thread::get_id();

            auto  threadIdValue =*(unsigned int *) &threadId;
            auto poolidx = threadIdValue % POOL_CNT;
            auto startidx=poolidx;
            int full_num = 0;
            RETRY:

            if(full_num>POOL_CNT){
                printf("run out of memory\n");
                exit(-1);
            }


            if(!alloc_mtx[poolidx].trylock()){
                if (blk_demand + cur_blk_[poolidx] > max_blk_[poolidx]) {
                    if(poolidx==startidx&&full_num<POOL_CNT){
                        full_num=0;
                    }
                    full_num++;
                }
                poolidx=(poolidx+1)%POOL_CNT;
                goto RETRY;
            }
            // case 1: not enough in the buffer
            if (blk_demand + cur_blk_[poolidx] > max_blk_[poolidx]) {
                if(poolidx==startidx&&full_num<POOL_CNT){
                    full_num=0;
                }
                full_num++;
                poolidx=(poolidx+1)%POOL_CNT;
                goto RETRY;
            }
            // case 2: current piece can not accommdate this allocation
            int piece_id = cur_blk_[poolidx] / piece_size_[poolidx];
            if ((cur_blk_[poolidx] % piece_size_[poolidx] + blk_demand) > piece_size_[poolidx]) {
                void *mem = buff_aligned_[poolidx][piece_id + 1]; // allocate from a new peice
                cur_blk_[poolidx] = piece_size_[poolidx] * (piece_id + 1) + blk_demand;
                meta_->metas[poolidx].cur_blk = cur_blk_[poolidx];
                Flush(&(meta_->metas[poolidx].cur_blk), 8);

                alloc_mtx[poolidx].unlock();
                return mem;
            }
                // case 3: current piece has enough space
            else {
                void *mem = buff_aligned_[poolidx][piece_id] + ALIGN_SIZE * (cur_blk_[poolidx] % piece_size_[poolidx]);

                cur_blk_[poolidx] = cur_blk_[poolidx] + blk_demand;
                meta_->metas[poolidx].cur_blk = cur_blk_[poolidx];
                Flush(&(meta_->metas[poolidx].cur_blk), 8);
                alloc_mtx[poolidx].unlock();
                return mem;
            }
        }

        static void free(void *addr) {
/*            for (int i = 0; i < PIECE_CNT; i++) {
                uint64_t offset = (uint64_t) addr - (uint64_t) buff_aligned_[i];
                if (offset > 0 && offset < piece_size_ * ALIGN_SIZE) {
                    // the addr is in this piece, do not reclaim it
                    return;
                }
            }*/

            // larger than 4KB, reclaim it
            uint64_t *header = (uint64_t *) ((uint64_t) addr - 8);
            uint64_t offset = *header;

            alloc_mtx[POOL_CNT].lock();
            auto oid_ptr = pmemobj_oid((void *) ((uint64_t) addr - offset));
            alloc_mtx[POOL_CNT].unlock();

            TOID(char) ptr_cpy;
            TOID_ASSIGN(ptr_cpy, oid_ptr);
            POBJ_FREE(&ptr_cpy);

        }

        /*
         *  Distinguish from virtual memory address and offset in the pool
         *  Each memory piece allocated from the pool has an in-pool offset, which remains unchanged
         *  until reclaiment. We cannot ensure that the pool file is mapped at the same position at
         *  any time, so it may locate at different virtual memory addresses next time.
         *
         *  So the rule is that, using virtual memory when doing normal operations like to DRAM
         *  space, using offset to store link relationship, for exmaple, next pointer in linklist
         */

        /*
         *  convert an offset to a virtual memory address
         */
        template<typename T>
        static inline T *absolute(T *pmem_offset) {
            if (pmem_offset == NULL)
                return NULL;
            return reinterpret_cast<T *>(reinterpret_cast<uint64_t>(pmem_offset) + reinterpret_cast<char *>(_pop));
        }

        /*
         *  convert a virtual memory address to an offset
         */
        template<typename T>
        static inline T *relative(T *pmem_direct) {
            if (pmem_direct == NULL)
                return NULL;
            return reinterpret_cast<T *>(reinterpret_cast<char *>(pmem_direct) - reinterpret_cast<char *>(_pop));
        }

    private:
        static void *mem_alloc(size_t nsize) {
            //statistic.allocnum++;
            PMEMoid tmp;
            alloc_mtx[POOL_CNT].lock();
            auto ret = pmemobj_zalloc(_pop, &tmp, nsize, TOID_TYPE_NUM(char));
            if (ret) {
                std::cout << "Fail logging: " << ret << "; Size = " << nsize << std::endl;
                LOG_FATAL("Allocate: Allocation Error in PMEMoid 11");
            }
            void *mem = pmemobj_direct(tmp);
            alloc_mtx[POOL_CNT].unlock();
            assert(mem != nullptr);
            return mem;
        }
    };

}


#endif //PLTREE_ALLOCATOR_H
