//
// Created by Cshuang on 2020/10/10.
//

#ifndef PLTREE_TEST_SEGMENT_H
#define PLTREE_TEST_SEGMENT_H

#include <cmath>
#include <cstdint>
#include <cmath>
#include <vector>
#include <cstring>
#include "common.h"
#include "builder.h"
#include <iostream>
#include "allocator.h"
#include <functional>

namespace plt {


    using Alloc = my_alloc::MemPoolAllocator;

    enum OpStatus {
        NotFound = 0,
        Find = 1,
        Deleted = 2,
        Inserted = 3,
        Updated = 4,
        Failed = 5,
        SMO = 6,
        lock = 7,
        Retry = 8
    };

    template<typename KeyType, typename ValueType, typename Alloc=my_alloc::MemPoolAllocator>//PMAllocator>
    class Segment;

    template<typename KeyType, typename ValueType, typename Alloc=my_alloc::MemPoolAllocator>
    class Bucket {
        typedef struct KV<KeyType, ValueType> KV;
    public:
        Bucket() {
            //for(int i=0;i<BUCKET_LEN)
        }

        int insert(KeyType key, ValueType value) {
            for (int i = 0; i < BUCKET_LEN; ++i) {
                if (kvs[i].first == 0) {
                    kvs[i].second = value;
                    kvs[i].first = key;
                    Alloc::Persist(&kvs[i], sizeof(KV));
                    return 2;
                }
            }
            return 3;
            //persist
        }

        bool update(int pos, KeyType key, ValueType value) {
            if (key == kvs[pos].first) {
                kvs[pos].second = value;
                Alloc::Persist(&kvs[pos].second, sizeof(ValueType));
                return true;
            }
            return false;
        }

        int find(KeyType key, ValueType &value) {
            for (int i = 0; i < BUCKET_LEN; ++i) {
                if (key == kvs[i].first) {
                    value = kvs[i].second;
                    return 1;
                }
            }
            return 0;
        }

        int erase(int pos, KeyType key) {
            if (key == kvs[pos].first) {
                kvs[pos].first = FREE;
                Alloc::Persist(&kvs[pos].first, sizeof(KeyType));
                return true;
            }
            return false;
        }

        // private:
        KV kvs[BUCKET_LEN];
    };

    template<typename KeyType, typename ValueType>
    class Entry {
    public:
        Entry() {
        }

        Bucket<KeyType, ValueType> buckets[BUCKET_PER_ENTRY];
    };

    template<typename KeyType>
    class EntryMeta {
    public:
        uint32_t epoch;
        uint32_t size;
    };

    template<typename KeyType, typename ValueType, typename Alloc=my_alloc::MemPoolAllocator>
    class Layer {
        typedef struct KV<KeyType, ValueType> KV;
        typedef Layer<KeyType, ValueType, Alloc> DirType;
        typedef Entry<KeyType, ValueType> EntryType;
        typedef Bucket<KeyType, ValueType, Alloc> BucketType;
        typedef EntryMeta<KeyType> EMetaType;

    private:
    public:
        EntryType *entrys;
        DirType *nextlevel;
        EntryFp *FPEntrys;
        size_t level;
        EMetaType *metas;

        Layer(const size_t &level, uint64_t &fpPtr, EntryFp *lastFPEntrys, innerSlot<KeyType> *_slot) : level(level),
                                                                                                        nextlevel(
                                                                                                                nullptr) {
            //// my_alloc::statistic.maxlevel = std::max(my_alloc::statistic.maxlevel, level);

            size_t size = 1 << (INIT_ENTRY_BIT + FANOUT * level);
            entrys = reinterpret_cast<EntryType *>(Alloc::PMZAlloc(size * sizeof(EntryType)));
            FPEntrys = NULL;
            my_alloc::allocMem(FPEntrys, (int) size);
            fpPtr = reinterpret_cast<uint64_t>(FPEntrys);
            for (int i = 0; i < size; ++i) {
                FPEntrys[i].entryPtr = reinterpret_cast<uintptr_t>(&entrys[i]);
            }
            if (level > 0 && lastFPEntrys != NULL) {
                int lastSize = 1 << (INIT_ENTRY_BIT + FANOUT * (level - 1));
                for (int i = 0; i < lastSize; ++i) {
                    lastFPEntrys[i].nextLevel = FPEntrys;
                }
            }
            metas = reinterpret_cast<EMetaType *>(Alloc::PMZAlloc(size * sizeof(EMetaType)));///
            Alloc::Persist(&entrys, sizeof(DirType));//persist
            _slot->addLevel();
            Alloc::Persist(&_slot->levelnum, sizeof(uint16_t));
        }

        void *operator new(size_t size) {
            void *p = Alloc::PMZAlloc(size);
            return p;
        }

        static uint32_t reverseBitsInRange(uint32_t num, int start, int end) {
            uint32_t mask = (((1 << (end - start + 1)) - 1) << start);  // 创建一个掩码来选取指定范围的位
            uint32_t bitsToReverse = (num & mask);  // 提取需要被反转的位
            bitsToReverse >>= start;  // 将提取的位右移至最低位
            uint32_t reversedBits = 0;
            int i = end - start + 1;
            while (i != 0) {
                reversedBits <<= 1;
                reversedBits |= (bitsToReverse & 1);
                bitsToReverse >>= 1;
                i--;
            }
            return reversedBits;
        }

        KV *getKVAt(int entry, int bucket, int idx) {
            return &entrys[entry].buckets[bucket].kvs[idx];
        }

        void getAllKVs(int &deleted, std::map<KeyType, ValueType> &result) {
            int entryNum = 1 << (INIT_ENTRY_BIT + FANOUT * level);
            for (int i = entryNum - 1; i >= 0; --i) {
                if (metas[i].size == 0)
                    continue;
                int n = metas[i].size % BUCKET_LEN == 0 ? metas[i].size / BUCKET_LEN - 1 : metas[i].size / BUCKET_LEN;
                for (int j = n; j >= 0; --j) {
                    int k = (j < n ? BUCKET_LEN - 1 : metas[i].size % BUCKET_LEN - 1);
                    for (; k >= 0; --k) {
                        KV *kv_ = getKVAt(i, j, k);
                        if (entrys[i].buckets[j].kvs[k].first != FREE) {
                            result[entrys[i].buckets[j].kvs[k].first] = entrys[i].buckets[j].kvs[k].second;
                            const auto [ret, success] = result.insert(
                                    std::make_pair<KeyType, ValueType>(std::move(kv_->first), std::move(kv_->second)));
                            if (success && kv_->second == DELETED) {
                                ++deleted;
                            }
                        }
                    }
                }
            }
        }

        void rangeScan(KeyType key, int num, int &deleted, std::map<KeyType, ValueType> &results) {
            int result_size = results.size();
            uint64_t upper_bound = UINT64_MAX;
            if (result_size > 0) {
                auto it = prev(results.end());
                while (it != results.begin() && it->second == DELETED) {
                    --it;
                }
                if (it->second != DELETED)
                    upper_bound = it->first;
            }
            int entryNum = 1 << (INIT_ENTRY_BIT + FANOUT * level);
            for (int i = entryNum - 1; i >= 0; --i) {
                if (metas[i].size == 0)
                    continue;
                int n = metas[i].size % BUCKET_LEN == 0 ? metas[i].size / BUCKET_LEN - 1 : metas[i].size / BUCKET_LEN;
                for (int j = n; j >= 0; --j) {
                    int k = (j < n ? BUCKET_LEN - 1 : metas[i].size % BUCKET_LEN - 1);
                    for (; k >= 0; --k) {
                        KV *kv_ = getKVAt(i, j, k);
                        if (kv_->first != FREE && kv_->first >= key &&
                            (result_size - deleted < num || kv_->first < upper_bound)) {
                            const auto [ret, success] = results.insert(
                                    std::make_pair<KeyType, ValueType>(std::move(kv_->first), std::move(kv_->second)));
                            if (success) {
                                ++result_size;
                                if (kv_->second == DELETED)
                                    ++deleted;
                            }
                            if (result_size - deleted > num) {
                                // Delete the largest element
                                auto it = prev(results.end());
                                while (it != results.begin() && it->second == DELETED) {
                                    --it;
                                }
                                if (it->second != DELETED) {
                                    if (it != results.begin()) {
                                        auto prev_it = prev(it);
                                        results.erase(it);
                                        while (prev_it != results.begin() && prev_it->second == DELETED) {
                                            --prev_it;
                                        }
                                        if (prev_it->second != DELETED)
                                            upper_bound = prev_it->first;
                                        --result_size;
                                    } else {
                                        results.erase(it);
                                        --result_size;
                                        upper_bound = UINT64_MAX;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        static size_t hash_key(const KeyType &key) {
            size_t hash = std::_Hash_bytes(&key, sizeof(KeyType), HASH_SEED);
            return hash;
        }

        void sortInEntry(int entryIdx, std::vector<KV> &kvs) {
            for (int i = 0; i < BUCKET_PER_ENTRY; ++i) {
                for (int j = 0; j < BUCKET_LEN; ++j) {
                    kvs.push_back(entrys[entryIdx].buckets[i].kvs[j]);
                }
            }
            std::sort(kvs.begin(), kvs.end(), [](KV &a, KV &b) {
                return a.first < b.first;
            });
        }

        static uint32_t rehashEntryIdx(uint32_t entryIdx, int targetlevel) {
            if (!ENABLE_LBM)
                return entryIdx >> (INIT_ENTRY_BIT + FANOUT * (targetlevel - 1));
            else {
                int i = 1 << FANOUT - 1;
                return entryIdx & i;
            }
        }

        static uint32_t getEntryIdx(int targetlevel, KeyType key) {
            uint32_t entryIdx;
            uint32_t size = 1 << (INIT_ENTRY_BIT + FANOUT * targetlevel);
            size_t key_hash = hash_key(key);
            entryIdx = key_hash & (size - 1);
            if (ENABLE_LBM) {
                // We hash-partition: Calculate the hash and clamp to table size
                entryIdx = reverseBitsInRange(entryIdx, 0, INIT_ENTRY_BIT + FANOUT * targetlevel - 1);
            }
            return entryIdx;
        }

        void mergeToNextLevel(int entryIdx, innerSlot<KeyType> *_slot) {
            if (nextlevel == nullptr) {
                uint64_t fpPtr = 0;
                nextlevel = new DirType(level + 1, fpPtr, FPEntrys, _slot);
                Alloc::Persist(&nextlevel, sizeof(DirType *));
            }
            // EntryType* entry=&entrys[entryIdx];
            KV *kvs[BUCKET_PER_ENTRY * (BUCKET_LEN) * (1 << FANOUT)] = {nullptr};
            uint32_t rehashed = 0;
            uint32_t bucketIdx = 0;
            uint32_t sizes[1 << FANOUT] = {0};
            while (rehashed < metas[entryIdx].size && bucketIdx < BUCKET_PER_ENTRY) {
                int toRehash = std::min((uint32_t) BUCKET_LEN, metas[entryIdx].size - rehashed);
                rehash(nextlevel->level, entryIdx, bucketIdx, toRehash, kvs, sizes);
                rehashed += toRehash;
                ++bucketIdx;
            }
            bulkInsertToNextLevel(kvs, sizes, level + 1, _slot);
            //clear the entry
            FPEntrys[entryIdx].clear();
            //FPEntrys[entryIdx].size = 0;
            metas[entryIdx].size = 0;
            Alloc::Persist(&metas[entryIdx].size, sizeof(uint32_t));
        }

        void rehash(uint32_t targetLevel, uint32_t entryIdx, uint32_t bucketIdx, uint32_t rehash, KV **kvs,
                    uint32_t *sizes) {
            uint32_t kvnum = (BUCKET_LEN) * BUCKET_PER_ENTRY;//(1<<FANOUT);
            BucketType *bucket = &entrys[entryIdx].buckets[bucketIdx];
            for (uint32_t i = 0; i < rehash; ++i) {
                uint32_t newEntryIdx = getEntryIdx(targetLevel, bucket->kvs[i].first);
                newEntryIdx = rehashEntryIdx(newEntryIdx,
                                             targetLevel);//newEntryIdx >> (INIT_ENTRY_BIT + FANOUT * (targetLevel - 1));
                auto first = kvs + (newEntryIdx * kvnum);
                auto last = kvs + (newEntryIdx * kvnum + sizes[newEntryIdx]);
                auto pos = std::lower_bound(first, last, &bucket->kvs[i], [](KV *s1, KV *s2) {
                    return s1->first < s2->first;
                });
                if (pos != last) {
                    auto offset = pos - kvs;
                    if (bucket->kvs[i].first == kvs[offset]->first)
                        kvs[offset]->second = bucket->kvs[i].second;
                    continue;
                }
                kvs[newEntryIdx * kvnum + sizes[newEntryIdx]] = &bucket->kvs[i];
                ++sizes[newEntryIdx];
            }
        }

        void bulkInsertToNextLevel(KV **kvs, uint32_t *sizes,
                                   size_t targetlevel, innerSlot<KeyType> *_slot) {
            uint32_t kvnum = (BUCKET_LEN) * BUCKET_PER_ENTRY;
            uint32_t entrySize = 1 << (FANOUT);
            for (uint32_t i = 0; i < entrySize; ++i) {
                if (sizes[i] == 0) {
                    continue;
                }
                uint32_t start = i * kvnum;
                auto kv = kvs[start];
                uint32_t entryIdx = getEntryIdx(nextlevel->level, kv->first);
                EMetaType *meta = nextlevel->metas + entryIdx;
                uint32_t size_ = meta->size;
                if (size_ + sizes[i] > BUCKET_PER_ENTRY * BUCKET_LEN) {
                    nextlevel->mergeToNextLevel(entryIdx, _slot);
                }
                uint32_t retNum = bulkInsert(nextlevel, entryIdx, sizes[i], start, kvs, targetlevel);
                assert (retNum == sizes[i]);
            }
        }

        int getFreeBucketIdx(uint32_t size) {
            int freeBucketIdx = size / BUCKET_LEN;
            return freeBucketIdx < BUCKET_PER_ENTRY ? freeBucketIdx : -1;
        }

        uint32_t getSizeOfBucket(uint32_t size, uint32_t bucketIdx) {
            if (size - bucketIdx * BUCKET_LEN > 0 && size - bucketIdx * BUCKET_LEN < BUCKET_LEN) {
                return size - bucketIdx * BUCKET_LEN;
            } else {
                return 0;
            }
        }

        int
        bulkInsert(DirType *dir, uint32_t nextLevelEntryIdx, uint32_t size, uint32_t start, KV **kvs,
                   size_t targetlevel) {
            EntryType *entry = dir->entrys + nextLevelEntryIdx;
            EMetaType *meta = dir->metas + nextLevelEntryIdx;
            int n = getFreeBucketIdx(meta->size);
            if (n == -1) {
                // All buckets are already full
                return 0;
            }
            uint32_t elemsInserted = 0;
            // EntryFp * fp_entrys = fpPtr);
            while (n < BUCKET_PER_ENTRY && elemsInserted < size) {
                BucketType *bucket = entry->buckets + n;
                uint32_t bucketSize = getSizeOfBucket(meta->size, n);
                auto fp_bucket = dir->FPEntrys[nextLevelEntryIdx].fps + n * BUCKET_LEN;
                uint32_t elemsToInsert = std::min(BUCKET_LEN - bucketSize, size - elemsInserted);
                if (elemsToInsert > 0) {
                    //  insert_into_filter(keys + elems_inserted, elems_to_insert, level, directory_entry_idx, n);
                    for (uint32_t i = 0; i < elemsToInsert; ++i) {
                        uint32_t offset = start + elemsInserted + i;
                        _mm_stream_si64((long long *) (&bucket->kvs[bucketSize + i].first),
                                        (long long) kvs[offset]->first);
                        _mm_stream_si64((long long *) (&bucket->kvs[bucketSize + i].second),
                                        (long long) kvs[offset]->second);
                    }
                    for (uint32_t i = 0; i < elemsToInsert; ++i) {
                        uint32_t offset = start + elemsInserted + i;
                        unsigned char fp = hashcode1B<KeyType>(kvs[offset]->first);
                        fp_bucket[bucketSize + i] = fp;
                    }
                    elemsInserted += elemsToInsert;
                }
                ++n;
            }
            //meta->epoch = epoch;
            meta->size += elemsInserted;
            //dir->FPEntrys[nextLevelEntryIdx].size += elemsInserted;
            Alloc::Persist(&meta->size, sizeof(uint32_t));
            return elemsInserted;
        }


    };

    template<typename KeyType, typename ValueType, typename Alloc>//PMAllocator>
    class Segment {
        typedef struct KV<KeyType, ValueType> KV;
        typedef Layer<KeyType, ValueType> DirType;
        typedef Entry<KeyType, ValueType> EntryType;
        typedef Bucket<KeyType, ValueType, Alloc> BucketType;
        typedef std::function<OpStatus(const KeyType &, ValueType &, const int &, blockMeta *,
                                       innerSlot<KeyType> *)> funtype;
    public:
        Segment() : num_array_keys_(0), slope_(0.0), num_buffers_keys_(0), pre_(nullptr),
                    next_(nullptr) {

        }

        ~Segment() {

        }

        /*     void initFunc(){
                 functions[0] =[this](const KeyType &key, ValueType &value, const int &idx, blockMeta *bm, innerSlot<KeyType> *_slot) -> OpStatus {
                     uint32_t old_version = __atomic_load_n(&blocks[idx].head.version_lock, __ATOMIC_ACQUIRE);
                     uint8_t pos = __atomic_load_n(&blocks[idx].head.num, __ATOMIC_ACQUIRE);
                     uint64_t buffer = blocks[idx].head.getBufferPtr();
                     if (old_version & lockSet) {
                         //sched_yield();
                         return OpStatus::lock;
                     }
                     OpStatus ret = NotFound;
                     int kvPerBlock = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
                     //uint8_t fp = hashcode1B<KeyType>(key);
                     for (int i = 0; i <kvPerBlock; ++i) {
                         if(blocks[idx].kv[i].first==FREE){
                             break;
                         }
                         if (key_equal(blocks[idx].kv[i].first ,key)) {
                             value = blocks[idx].kv[i].second;
                             if (value == DELETED){
                                 ret=OpStatus::Deleted;
                             }
                             else{
                                 ret=OpStatus::Find;
                             }
                             break;
                         }
                     }
                     if(ret==OpStatus::Deleted){
                         if (blocks[idx].head.check_version_change(old_version) || blocks[idx].head.check_size_change(pos))
                             return OpStatus::Retry;
                         return OpStatus::NotFound;
                     }
                     else if(ret==OpStatus::Find){
                         if (blocks[idx].head.check_version_change(old_version) || blocks[idx].head.check_size_change(pos))
                             return OpStatus::Retry;
                         return OpStatus::Find;
                     }

                     auto ret1=false;

                     if (unlikely(buffer != 0)) {
                         auto meta=bm+idx;
                         ret1 = searchInDir(key, value, reinterpret_cast<EntryFp *>(meta->getFPtr()));
                         if (meta->check_version_change(old_version) || meta->check_size_change(pos))
                             return OpStatus::Retry;
                     }
                     return ret1?OpStatus::Find:OpStatus::NotFound;
                 };
                 functions[1] =[this](const KeyType &key, ValueType &value, const int &idx, blockMeta *bm, innerSlot<KeyType> *_slot) -> OpStatus {
                     auto meta=bm+idx;
                     if(unlikely(meta->v!=_slot->v)){
                         if(!meta->try_get_lock()){
                             return OpStatus::Retry;
                         }
                         rebuildBlockMeta(idx,_slot->v,meta);
                         meta->release_lock_without_change();
                     }
                     uint32_t old_version = __atomic_load_n(&meta->version_lock, __ATOMIC_ACQUIRE);
                     uint8_t pos = __atomic_load_n(&meta->num, __ATOMIC_ACQUIRE);
                     uint64_t buffer = meta->getLayerPtr();

                     if (old_version & lockSet) {
                         //sched_yield();
                         return OpStatus::lock;
                     }
                     OpStatus ret = NotFound;
                     int kvPerBlock = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
                     uint8_t fp = hashcode1B<KeyType>(key);
                     for (int i = 0; i <kvPerBlock; ++i) {
                         if (meta->fps[i] == fp) {
                             if(key_equal(blocks[idx].kv[i].first ,key)){
                                 value = blocks[idx].kv[i].second;
                                 if (value == DELETED)
                                     ret=OpStatus::Deleted;
                                 else
                                     ret=OpStatus::Find;
                                 break;
                             }
                         }
                     }

                     if(ret==OpStatus::Deleted){
                         if (meta->check_version_change(old_version) || meta->check_size_change(pos))
                             return OpStatus::Retry;
                         return OpStatus::NotFound;
                     }
                     else if(ret==OpStatus::Find) {
                         if (meta->check_version_change(old_version) || meta->check_size_change(pos))
                             return OpStatus::Retry;
                         return OpStatus::Find;
                     }

                     auto ret1=false;
                     if (likely(buffer != 0)) {
                         ret1 = searchInDir(key, value, reinterpret_cast<EntryFp *>(meta->getFPtr()));
                         if (meta->check_version_change(old_version) || meta->check_size_change(pos))
                             return OpStatus::Retry;
                     }
                     return ret1?OpStatus::Find:OpStatus::NotFound;
                 };
             }*/
        void rebuildBlockMeta(int idx, uint8_t version, blockMeta *bm) {
            uint8_t kvPerBlock = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
            uint8_t num = 0;
            auto ptr = blocks[idx].head.ptr;
            bm->setLayerPtr(ptr);
            bm->v = version;
            for (uint8_t i = 0; i < kvPerBlock; ++i) {
                if (blocks[idx].kv[i].first != FREE) {
                    uint8_t fp = hashcode1B<KeyType>(blocks[idx].kv[i].first);
                    bm->fps[i] = fp;
                    bm->num++;
                    //blocks[idx].head.num++;
                } else {
                    break;
                }
            }
            EntryFp *lastFPEntrys = NULL;
            auto layer = reinterpret_cast< DirType *> (ptr);
            int level = 0;
            while (layer) {
                auto entrys = layer->entrys;
                auto metas = layer->metas;
                size_t size = 1 << (INIT_ENTRY_BIT + FANOUT * level);
                my_alloc::allocMem(layer->FPEntrys, (int) size);

                for (int i = 0; i < size; ++i) {
                    layer->FPEntrys[i].entryPtr = reinterpret_cast<uintptr_t>(&entrys[i]);
                    int count = 0;
                    bool breakFlag = false;
                    for (int j = 0; j < BUCKET_PER_ENTRY; ++j) {
                        auto bucket = entrys[i].buckets + j;
                        auto fp_bucket = layer->FPEntrys[i].fps + j * BUCKET_LEN;
                        for (int k = 0; k < BUCKET_LEN; ++k) {
                            if (bucket->kvs[k].first != FREE) {
                                uint8_t fp = hashcode1B<KeyType>(bucket->kvs[k].first);
                                fp_bucket[k] = fp;
                                ++count;
                            }
                            if (count == metas[i].size) {
                                breakFlag = true;
                                break;
                            }
                        }
                        if (breakFlag) {
                            break;
                        }
                    }
                }
                if (level == 0) {
                    bm->setFPtr(reinterpret_cast<uint64_t>(layer->FPEntrys));
                } else if (level > 0 && lastFPEntrys != NULL) {
                    int lastSize = 1 << (INIT_ENTRY_BIT + FANOUT * (level - 1));
                    for (int i = 0; i < lastSize; ++i) {
                        lastFPEntrys[i].nextLevel = layer->FPEntrys;
                    }

                }
                lastFPEntrys = layer->FPEntrys;
                layer = layer->nextlevel;
            }
        }

        OpStatus
        upsert(const KeyType &key, const ValueType &value, const int &idx, blockMeta *bm, innerSlot<KeyType> *_slot) {

            uint8_t kvPerBlock = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
            if (idx == -1 || _slot == NULL)
                return OpStatus::Failed;
            if (!_slot->check_read_lock())
                return OpStatus::lock;
            if (_slot->levelnum / _slot->size > THRESHOLD)
                return OpStatus::SMO;
            auto meta = bm + idx;
            if (!_slot->check_write_lock())
                return OpStatus::lock;
            meta->get_lock();
            if (meta->v != _slot->v) {
                rebuildBlockMeta(idx, _slot->v, meta);
                // blocks[idx].head.v=_slot->v;
                //blocks[idx].head.release_lock_without_change();
            }
            // blocks[idx].head.get_lock();
            RETRY:
            uint8_t pos = meta->num;
            // uint8_t pos2 = blocks[idx].head.num;
            if (pos >= kvPerBlock) {
                migrate(idx, bm, _slot);
                // blocks[idx].head.clearNum();
                goto RETRY;
            }
            uint8_t fp = hashcode1B<KeyType>(key);
            for (uint8_t i = 0; i < pos; ++i) {
                if (meta->fps[i] == fp) {
                    if (key_equal(blocks[idx].kv[i].first, key)) {
                        blocks[idx].kv[i].second = value;
                        Alloc::Persist(&blocks[idx].kv[i].second, sizeof(ValueType));
                        meta->release_lock();
                        // blocks[idx].head.release_lock();
                        return OpStatus::Updated;
                    }
                }
            }
            blocks[idx].kv[pos].second = value;
            blocks[idx].kv[pos].first = key;
            meta->fps[pos] = fp;
            Alloc::Persist(&blocks[idx].kv[pos], sizeof(KV));
            meta->addNum();
            //blocks[idx].head.addNum();
            meta->release_lock();
            // blocks[idx].head.release_lock();
            return OpStatus::Inserted;

        }

        inline bool key_equal(const KeyType &a, const KeyType &b) const {
            return !(a < b) && !(b < a);
        }
        OpStatus
        upsert(const KeyType &key, const ValueType &value,  innerSlot<KeyType> *_slot) {
            uint8_t kvPerBlock = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
            auto bm=reinterpret_cast<blockMeta * >(_slot->getBufferPtr());
            auto block=(key-_slot->minkey) * _slot->slope /kvPerBlock;
            int idx=block<0?0:(block>=_slot->size?_slot->size-1:block);
            if (idx == -1 || _slot == NULL)
                return OpStatus::Failed;
            if (!_slot->check_read_lock())
                return OpStatus::lock;
            if (_slot->levelnum / _slot->size > THRESHOLD)
                return OpStatus::SMO;
            auto meta = bm + idx;
            if (!_slot->check_write_lock())
                return OpStatus::lock;
            meta->get_lock();
            if (meta->v != _slot->v) {
                rebuildBlockMeta(idx, _slot->v, meta);
                // blocks[idx].head.v=_slot->v;
                //blocks[idx].head.release_lock_without_change();
            }
            // blocks[idx].head.get_lock();
            RETRY:
            uint8_t pos = meta->num;
            // uint8_t pos2 = blocks[idx].head.num;
            if (pos >= kvPerBlock) {
                migrate(idx, bm, _slot);
                // blocks[idx].head.clearNum();
                goto RETRY;
            }
            uint8_t fp = hashcode1B<KeyType>(key);
            for (uint8_t i = 0; i < pos; ++i) {
                if (meta->fps[i] == fp) {
                    if (key_equal(blocks[idx].kv[i].first, key)) {
                        blocks[idx].kv[i].second = value;
                        Alloc::Persist(&blocks[idx].kv[i].second, sizeof(ValueType));
                        meta->release_lock();
                        // blocks[idx].head.release_lock();
                        return OpStatus::Updated;
                    }
                }
            }
            blocks[idx].kv[pos].second = value;
            blocks[idx].kv[pos].first = key;
            meta->fps[pos] = fp;
            Alloc::Persist(&blocks[idx].kv[pos], sizeof(KV));
            meta->addNum();
            //blocks[idx].head.addNum();
            meta->release_lock();
            // blocks[idx].head.release_lock();
            return OpStatus::Inserted;

        }

        OpStatus find(const KeyType &key, ValueType &value, innerSlot<KeyType> *_slot) {
            auto bm=reinterpret_cast<blockMeta * >(_slot->getBufferPtr());
            auto block=(key-_slot->minkey) * _slot->slope /15;
            int idx=block<0?0:(block>=_slot->size?_slot->size-1:block);
            if (idx == -1)
                return OpStatus::Failed;
            if (!_slot->check_read_lock()) {
                return OpStatus::lock;
            }
           // auto meta = bm + idx;
//            if (__glibc_unlikely(meta->v != _slot->v)) {
//                if (!meta->try_get_lock()) {
//                    return OpStatus::Retry;
//                }
//                rebuildBlockMeta(idx, _slot->v, meta);
//                meta->release_lock_without_change();
//            }
auto head=&blocks[idx].head;
            uint32_t old_version = __atomic_load_n(&head->version_lock, __ATOMIC_ACQUIRE);
            uint8_t pos = __atomic_load_n(&head->num, __ATOMIC_ACQUIRE);
            if (old_version & lockSet) {
                //sched_yield();
                return OpStatus::lock;
            }
            //OpStatus ret = NotFound;
            auto kvPerBlock = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
            //uint8_t fp = hashcode1B<KeyType>(key);
            for (int i = 0; i < kvPerBlock; ++i) {
                if (key_equal(blocks[idx].kv[i].first, key)) { //meta->fps[i] == fp &&
                    value = blocks[idx].kv[i].second;
                    if (head->check_version_change(old_version) || head->check_size_change(pos))
                        return OpStatus::Retry;
                    return OpStatus::Find;
                }
            }
            auto ret1 = false;
            uint64_t buffer = head->ptr;
            if (__glibc_likely(buffer!=0)) {
                //
                auto meta = bm + idx;
                ret1 = searchInDir(key, value, reinterpret_cast<EntryFp *>(meta->getFPtr()));
                if (head->check_version_change(old_version) || head->check_size_change(pos))
                    return OpStatus::Retry;
            }
            return ret1 ? OpStatus::Find : OpStatus::NotFound;
        }
        OpStatus find(const KeyType &key, ValueType &value, const int &idx, blockMeta *bm, innerSlot<KeyType> *_slot) {
            if (idx == -1)
                return OpStatus::Failed;
            if (!_slot->check_read_lock()) {
                return OpStatus::lock;
            }
            auto meta = bm + idx;
            if (__glibc_unlikely(meta->v != _slot->v)) {
                if (!meta->try_get_lock()) {
                    return OpStatus::Retry;
                }
                rebuildBlockMeta(idx, _slot->v, meta);
                meta->release_lock_without_change();
            }
            uint32_t old_version = __atomic_load_n(&meta->version_lock, __ATOMIC_ACQUIRE);
            uint8_t pos = __atomic_load_n(&meta->num, __ATOMIC_ACQUIRE);
            if (old_version & lockSet) {
                //sched_yield();
                return OpStatus::lock;
            }
            //OpStatus ret = NotFound;
            auto kvPerBlock = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
            uint8_t fp = hashcode1B<KeyType>(key);
            for (int i = 0; i < kvPerBlock; ++i) {
                if (meta->fps[i] == fp && key_equal(blocks[idx].kv[i].first, key)) {
                    value = blocks[idx].kv[i].second;
                    if (meta->check_version_change(old_version) || meta->check_size_change(pos))
                        return OpStatus::Retry;
                    return OpStatus::Find;
                }
            }
            auto ret1 = false;
            if (__glibc_likely(meta->hasBuffer())) {
                //uint64_t buffer = meta->getLayerPtr();
                ret1 = searchInDir(key, value, reinterpret_cast<EntryFp *>(meta->getFPtr()));
                if (meta->check_version_change(old_version) || meta->check_size_change(pos))
                    return OpStatus::Retry;
            }
            return ret1 ? OpStatus::Find : OpStatus::NotFound;
        }
/*         {
             if(idx==-1)
                 return OpStatus::Failed;
             if(!_slot->check_read_lock()){
                 return OpStatus::lock;
             }
             if (_slot->levelnum/_slot->size<=THRESHOLD1){
                 auto meta=bm+idx;
                 if(blocks[idx].head.v!=_slot->v){
                     if(!meta->try_get_lock()){
                         return OpStatus::Retry;
                     }
                     rebuildBlockMeta(idx,_slot->v,meta);
                     blocks[idx].head.v=_slot->v;
                     blocks[idx].head.release_lock_without_change();
                     meta->release_lock_without_change();
                 }
                 uint32_t old_version = __atomic_load_n(&blocks[idx].head.version_lock, __ATOMIC_ACQUIRE);
                 uint8_t pos = __atomic_load_n(&blocks[idx].head.num, __ATOMIC_ACQUIRE);*//* *//*
                 *//*  *//*
                 if (old_version & lockSet) {
                     //sched_yield();
                     return OpStatus::lock;
                 }
                 OpStatus ret = NotFound;
                 int kvPerBlock = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
                 //uint8_t fp = hashcode1B<KeyType>(key);
                 for (int i = 0; i <kvPerBlock; ++i) {
                     if (key_equal(blocks[idx].kv[i].first ,key)) {
                         value = blocks[idx].kv[i].second;
                         if (value == DELETED)
                             ret=OpStatus::Deleted;
                         else
                             ret=OpStatus::Find;
                     }
                 }
                 if(ret==OpStatus::Deleted){
                     if (blocks[idx].head.check_version_change(old_version) || blocks[idx].head.check_size_change(pos))
                         return OpStatus::Retry;
                     return OpStatus::NotFound;
                 }
                 else if(ret==OpStatus::Find){
                     if (blocks[idx].head.check_version_change(old_version) || blocks[idx].head.check_size_change(pos))
                         return OpStatus::Retry;
                     return OpStatus::Find;
                 }

                 auto ret1=false;
                 uint64_t buffer = blocks[idx].head.getBufferPtr();
                 if (buffer != 0) {
                     ret1 = searchInDir(key, value, reinterpret_cast<EntryFp *>(meta->getFPtr()));
                     if (meta->check_version_change(old_version) || meta->check_size_change(pos))
                         return OpStatus::Retry;
                 }
                 return ret1?OpStatus::Find:OpStatus::NotFound;
             }
             else
             {
                 auto meta=bm+idx;
                 if(meta->v!=_slot->v){
                     if(!meta->try_get_lock()){
                         return OpStatus::Retry;
                     }
                     rebuildBlockMeta(idx,_slot->v,meta);
                     meta->release_lock_without_change();
                 }

                 uint32_t old_version = __atomic_load_n(&meta->version_lock, __ATOMIC_ACQUIRE);
                 uint8_t pos = __atomic_load_n(&meta->num, __ATOMIC_ACQUIRE);*//* *//*
                 *//*  *//*
                 if (old_version & lockSet) {
                     //sched_yield();
                     return OpStatus::lock;
                 }
                 OpStatus ret = NotFound;
                 int kvPerBlock = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
                 uint8_t fp = hashcode1B<KeyType>(key);
                 for (int i = 0; i <kvPerBlock; ++i) {
                     if (meta->fps[i] == fp&&key_equal(blocks[idx].kv[i].first ,key)) {
                         value = blocks[idx].kv[i].second;
                         if (value == DELETED)
                             ret=OpStatus::Deleted;
                         else
                             ret=OpStatus::Find;
                     }
                 }
                 if(ret==OpStatus::Deleted){
                     if (meta->check_version_change(old_version) || meta->check_size_change(pos))
                         return OpStatus::Retry;
                     return OpStatus::NotFound;
                 }
                 else if(ret==OpStatus::Find){
                     if (meta->check_version_change(old_version) || meta->check_size_change(pos))
                         return OpStatus::Retry;
                     return OpStatus::Find;
                 }
                 auto ret1=false;
                 uint64_t buffer = meta->getLayerPtr();
                 if (buffer != 0) {
                     ret1 = searchInDir(key, value, reinterpret_cast<EntryFp *>(meta->getFPtr()));
                     if (meta->check_version_change(old_version) || meta->check_size_change(pos))
                         return OpStatus::Retry;
                 }
                 return ret1?OpStatus::Find:OpStatus::NotFound;
             }
        }*/

/*        inline bool find(KeyType key,ValueType &value,SegInfType * seginf){
            int idx = predict_block(key,seginf);
            auto blockHead= &(blocks+idx)->head;
            //prefetch(idx);
            RETRY:
            uint32_t old_version =
                    __atomic_load_n(&blockHead->version_lock, __ATOMIC_ACQUIRE);
            uint8_t pos=__atomic_load_n(&blockHead->num, __ATOMIC_ACQUIRE);
            if (old_version & lockSet) {
                sched_yield();
                goto RETRY;
            }
            if(pos!=0){
                auto ret=searchInBlock(key,value,idx);
                if(blockHead->check_version_change(old_version)||blockHead->check_size_change(pos))
                    goto RETRY;
                if(ret==OpStatus::Deleted)
                    return 0;
                else if(ret==OpStatus::Find)
                    return 1;
            }
            return searchInDir(key,value,idx);

        }*/
        inline OpStatus erase(const KeyType &key, const int &idx, blockMeta *bm, innerSlot<KeyType> *_slot) {
            return upsert(key, reinterpret_cast<ValueType>(DELETED), idx, bm, _slot);
            // return ret;
        }
        inline OpStatus erase(const KeyType &key,  innerSlot<KeyType> *_slot) {
            return upsert(key, reinterpret_cast<ValueType>(DELETED),  _slot);
            // return ret;
        }

        OpStatus searchInBlock(KeyType key, ValueType &value, int idx) {
            uint8_t kvPerBlock = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
            for (int i = kvPerBlock - 1; i >= 0; --i) {
                if (blocks[idx].kv[i].first == key) {
                    value = blocks[idx].kv[i].second;
                    if (value == DELETED)
                        return OpStatus::Deleted;
                    else
                        return OpStatus::Find;
                }
            }
            return OpStatus::NotFound;
        }

        bool findInLayer(KeyType key, ValueType &value, int n, EntryType *_entry) {
            ////my_alloc::statistic.bufferused++;
            int bucketIdx = n / BUCKET_LEN;
            BucketType *bucket = &_entry->buckets[bucketIdx];
            if (bucket->kvs[n % BUCKET_LEN].first == key) {
                value = bucket->kvs[n % BUCKET_LEN].second;
                return true;
            }
            return false;
        }


        bool searchInDir(KeyType key, ValueType &value, EntryFp *fp_entrys) {
            //int kvnum = BUCKET_LEN * BUCKET_PER_ENTRY;
            unsigned char fp = hashcode1B<KeyType>(key);
            __m128i target = _mm_set1_epi8(fp);
            int level = 0;
            while (fp_entrys != NULL) {
                int entryIdex = DirType::getEntryIdx(level, key);
                for (int i = BUCKET_PER_ENTRY - 1; i >= 0; i--) {
                    __m128i data = _mm_loadu_si128((__m128i *) &fp_entrys[entryIdex].fps[i * 16]);
                    __m128i compare = _mm_cmpeq_epi8(data, target);
                    int mask = _mm_movemask_epi8(compare);
                    // int pos2= __builtin_ctz(1);
                    while (mask != 0) {
                        int pos1 = __builtin_ctz(mask);
                        int position = i * 16 + pos1;
                        auto ret = findInLayer(key, value, position,
                                               reinterpret_cast<EntryType *>(fp_entrys[entryIdex].entryPtr));
                        if (ret) {
                            if (value != DELETED)
                                return true;
                            else
                                return false;
                        }
                        mask = mask & ~(1 << pos1);
                    }
                }
                level++;
                fp_entrys = fp_entrys[entryIdex].nextLevel;
            }
            return false;
        }

        void migrate(const int &idx, blockMeta *bm, innerSlot<KeyType> *_slot) {
            auto blockHead = &blocks[idx].head;
            auto dir = reinterpret_cast<DirType *> (blockHead->getBufferPtr());
            newDir(dir, idx, 0, bm, _slot);
            ///auto fpPtr=blockHead->getFpPtr();
            KV *kvs[1 * BUCKET_LEN * (1 << INIT_ENTRY_BIT)] = {
                    nullptr};//数组长度为分配到下一层level0的每个entry最多kv数总和=本层一个entry的kv数*fanout
            uint32_t sizes[1 << (INIT_ENTRY_BIT)] = {0}; //level0的entry数
            uint32_t num = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
            uint32_t kvnum = 1 * BLOCK_SIZE / sizeof(KV);//上一层的kv数/entry
            for (int i = 0; i < num; ++i) {
                uint32_t newEntryIdx = DirType::getEntryIdx(0, getKeyAt(idx, i));
                kvs[newEntryIdx * kvnum + sizes[newEntryIdx]] = blocks[idx].kv + i;
                sizes[newEntryIdx]++;
            }
            int entrySize = 1 << INIT_ENTRY_BIT;
            for (int i = 0; i < entrySize; ++i) {
                if (sizes[i] == 0) {
                    continue;
                }
                uint32_t start = i * kvnum;
                auto key = kvs[start]->first;
                uint32_t entryIdx = DirType::getEntryIdx(0, key);
                //EntryType* entry=dir->entrys+entryIdx;
                EntryMeta<KeyType> *meta = dir->metas + entryIdx;
                if (meta->size + sizes[i] > BUCKET_PER_ENTRY * BUCKET_LEN) {
                    dir->mergeToNextLevel(entryIdx, _slot);
                }
                int retNum = dir->bulkInsert(dir, entryIdx, sizes[i], start, kvs, 0);
                assert (retNum == sizes[i]);
            }
            // blockHead->clearNum();

            for (int i = 0; i < num; ++i) {
                delKeyUnsafeAt(idx, i);
            }
            Alloc::Persist(&blocks[idx], BLOCK_SIZE);
            bm[idx].clearNum();
            bm[idx].clear();
            //blockHead->epoch++;///
        }


        void updateResults(const KeyType &lower_bound, std::vector<KV> &kvs, const int &num, const int &deleted,
                           std::map<KeyType, ValueType> &results) {
            int result_size = results.size();
            uint64_t upper_bound = UINT64_MAX;
            if (result_size > 0) {
                auto it = prev(results.end());
                while (it != results.begin() && it->second == DELETED) {
                    --it;
                }
                if (it->second != DELETED)
                    upper_bound = it->first;
            }
            for (int i = 0; i < kvs.size(); ++i) {
                if (kvs[i].first >= lower_bound &&
                    (result_size - deleted < num || kvs[i].first < upper_bound)) {
                    const auto [it, success] = results.insert(
                            std::make_pair<KeyType, ValueType>(std::move(kvs[i].first), std::move(kvs[i].second)));
                    if (success) {
                        ++result_size;
                    }
                    if (result_size - deleted > num) {
                        // Delete the largest element
                        results.erase(prev(results.end()));
                        upper_bound = prev(results.end())->first;
                        --result_size;
                    }
                }
            }
        }

        int getAllKVs(std::map<KeyType, ValueType> &results, const uint32_t &blocknum) {
            int blockSize = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
            int deleted = 0;
            for (int idx = 0; idx < blocknum; ++idx) {
                for (int i = 0; i < blockSize; ++i) {
                    auto kv = getKVAt(idx, i);
                    if (kv.first != FREE) {
                        const auto [it, success] = results.insert(
                                std::make_pair<KeyType, ValueType>(std::move(kv.first), std::move(kv.second)));
                        (success && kv.second == DELETED) ? deleted++ : deleted;
                    }
                }
                if (blocks[idx].head.ptr != 0) {
                    auto ptr = reinterpret_cast<DirType *>(blocks[idx].head.ptr);
                    while (ptr != nullptr) {
                        ptr->getAllKVs(deleted, results);
                        ptr = ptr->nextlevel;
                    }
                }
            }
            return results.size();

        }

        int
        rangeScan(const KeyType &key, const int &num, int &deleted, std::map<KeyType, ValueType> &results, int idx,
                  const int &blocknum) {
            if (idx == -1)
                return 0;
            int ret = scanAtBlock(idx, key, num, deleted, results, blocknum);
            idx++;
            while (ret < num&&idx < blocknum  ) {
                ret = scanAtBlock(idx, key, num, deleted, results, blocknum);
                idx++;
            }
            return ret;
        }


        int
        scanAtBlock(const int &idx, const KeyType &key, const int &num, int &deleted,
                    std::map<KeyType, ValueType> &results, const int &blocknum) {
            std::vector<KV> kvs;
            //LOG("a");
            int blockSize = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
            if (blocknum > 0) {
                for (int i = 0; i < blockSize; ++i) {
                    auto kv = getKVAt(idx, i);
                    if (kv.first != FREE && kv.first >= key) {
                        kvs.push_back(kv);
                        kv.second == DELETED ? deleted++ : deleted;
                    }
                }
            }

            updateResults(key, kvs, num, deleted, results);
           // LOG("b");
            if (blocks[idx].head.ptr != 0) {
                auto ptr = reinterpret_cast<DirType *>(blocks[idx].head.ptr);
                while (ptr != nullptr) {
                    ptr->rangeScan(key, num, deleted, results);
                    ptr = ptr->nextlevel;
                }
            }
            //LOG("c");
            return results.size() - deleted;
        }

        int
        scanAtBegin(KeyType key, int num, int &deleted, std::map<KeyType, ValueType> &results, uint32_t blocknum) {
            int idx = 0;
            int ret = scanAtBlock(idx, key, num, deleted, results, blocknum);
            while (idx < blocknum && ret < num) {
                idx++;
                ret = scanAtBlock(idx, key, num, deleted, results, blocknum);
            }
            return ret;
        }

        inline OpStatus BulkloadKV(KeyType key, ValueType value, int idx) {
            uint8_t kvPerBlock = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
            for (int i = 0; i < kvPerBlock; --i) {
                if (blocks[idx].kv[i].first == FREE) {
                    blocks[idx].kv[i].second = value;
                    blocks[idx].kv[i].first = key;
                    return OpStatus::Inserted;
                }
            }
            return OpStatus::Failed;
        }

        inline int predict_block(KeyType key) const {
            if (key < minKey) {
                return 0;
            }
            int predicted_block = (key - minKey) * slope_ / ((BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV));
            if (predicted_block < 0)
                return 0;
            else if (predicted_block > blockNum - 1)
                return blockNum - 1;
            return predicted_block;
        }


        inline int predict_block(const KeyType &key, const uint32_t &block_number) const {

            int predicted_block = (key - minKey) * slope_ / ((BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV));
            if (predicted_block < 0)
                return 0;
            else if (predicted_block > block_number - 1)
                return block_number - 1;
            return predicted_block;
        }

        inline void
        newDir(DirType *&dir, const int &idx, const size_t &level, blockMeta *bm, innerSlot<KeyType> *_slot) {
            if (dir == nullptr) {
                uint64_t fpPtr = 0;
                dir = new DirType(level, fpPtr, NULL, _slot);
                addBufferAt(idx, (uint64_t) dir, fpPtr, bm);
            }
        }

        void freeMetas() {
            for (int i = 0; i < blockNum; ++i) {
                metas[i].freeFp();
            }
#ifndef LARGE_DRAM
            free(metas);
#endif
        }

        void init(uint32_t bnum, uint32_t keysNum) {
            blockNum = bnum;
            num_array_keys_ = keysNum;
            // initFunc();
        }

        void resize(const std::map<KeyType, ValueType> kvs, innerSlot<KeyType> *_slot, Segment<KeyType, ValueType> *pre,
                    Segment<KeyType, ValueType> *next) {
            int kvsPerBlock = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
            pre_ = pre;
            next_ = next;
            auto begin = kvs.begin();
            while (begin != kvs.end() && begin->second == DELETED) {
                begin++;
            }
            assert(begin != kvs.end());
            auto end = prev(kvs.end());
            while (end != kvs.begin() && end->second == DELETED) {
                --end;
            }
            assert(end != kvs.begin());
            slope_ = slope_ = num_array_keys_ * 1.0 / (end->first - begin->first);
            minKey = begin->first;
            slot = _slot;
            my_alloc::allocMem(metas, blockNum);
            slot->setBufferPtr(reinterpret_cast<uint64_t>(metas));
            while (begin != kvs.end()) {
                int idx = predict_block(begin->first, blockNum);
                auto block = blocks + idx;
                auto blockHead = &block->head;
                auto meta = metas + idx;
                uint8_t fp = hashcode1B(begin->first);
                RETRY:
                if (meta->num < kvsPerBlock) {
                    block->kv[meta->num].first = begin->first;
                    block->kv[meta->num].second = begin->second;
                    meta->fps[meta->num] = fp;
                    meta->num++;
                } else {
                    ////my_alloc::statistic.maxoverflow++;
                    migrate(idx, metas, _slot);
                    goto RETRY;
                }
                begin++;
            }

            return;
        }

        void bulk_loading(const std::vector<KV> &kvs, int start, innerSlot<KeyType> *_slot, uint8_t version) {

            if (num_array_keys_ == 0) {
                slope_ = 0;
                minKey = 0;
                slot = _slot;
                my_alloc::statistic.block0++;
                my_alloc::allocMem(metas, blockNum);
                for (int i = 0; i < blockNum; ++i) {
                    metas[i].v = version;
                    blocks[i].head.v = version;
                }
                slot->setBufferPtr(reinterpret_cast<uint64_t>(metas));
                /*  my_alloc::allocMem(metas, blockNum);
                *  for (int i = 0; i < blockNum; ++i) {
                       metas[i] = new std::vector<LayerMeta>();
                   }*/
            } else if (num_array_keys_ == 1) {
                slope_ = 0;
                minKey = kvs[start].first;
                slot = _slot;
                my_alloc::allocMem(metas, blockNum);
                for (int i = 0; i < blockNum; ++i) {
                    metas[i].v = version;
                    blocks[i].head.v = version;
                }
                slot->setBufferPtr(reinterpret_cast<uint64_t>(metas));
                auto block = blocks;
                auto blockHead = &block->head;
                auto meta = &metas[0];
                block->kv[meta->num] = kvs[start];
                meta->fps[meta->num] = hashcode1B(kvs[start].first);
                //blockHead->num++;
                meta->num++;
                // blockHead->num++;
                my_alloc::statistic.block15++;
/*                my_alloc::allocMem(metas, blockNum);
                for (int i = 0; i < blockNum; ++i) {
                    metas[i] = new std::vector<LayerMeta>();
                }*/

            } else if (num_array_keys_ >= 2) {
                int kvsPerBlock = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
                if (kvs[start + num_array_keys_ - 1].first == kvs[start].first) {
                    slope_ = 0;
                    minKey = kvs[start].first;
                } else {
                    slope_ = num_array_keys_ * 1.0 / (kvs[start + num_array_keys_ - 1].first - kvs[start].first);
                    minKey = kvs[start].first;
                }
                slot = _slot;
                my_alloc::allocMem(metas, blockNum);
                for (int i = 0; i < blockNum; ++i) {
                    metas[i].v = version;
                    blocks[i].head.v = version;
                }
                slot->setBufferPtr(reinterpret_cast<uint64_t>(metas));
                //minKey = kvs[start].first;

                int end = num_array_keys_ + start;

                /* double mean=0;
                 double st=0;
                 var_mean(kvs,start,num_array_keys_,mean,st);
                 int start_temp=start;
                 int end_temp=end;
                 bool flag=true;
                 for(int i=start;i<end;++i){
                     if (std::abs(kvs[i].first - mean) >  st) {
                         if(flag){
                             start_temp=i;
                         }else{
                             end_temp=i;
                             break;
                         }
                     }else
                         flag=false;
                 }
                 minKey=kvs[start_temp].first;
                 maxKey=kvs[end_temp].first;
                 slope_ = num_array_keys_*1.0 / (maxKey -minKey);*/
/*                auto it=std::lower_bound(kvs.begin()+start, kvs.begin()+end, KV{mean,0},[](const KV &kv1, const KV &kv2) {
                    return kv1.first < kv2.first;
                });
                int left=it-kvs.begin()-start;
                int right=num_array_keys_-left;
                if(left<right){
                    slope_=num_array_keys_*1.0 / (kvs[left*2+start].first - kvs[start].first);
                }else{
                    slope_=num_array_keys_*1.0 / (kvs[end-1].first - kvs[num_array_keys_-right*2+start].first);
                    minKey=kvs[num_array_keys_-right*2+start].first;
                }*/
                slope_ = slope_ / INIT_RATIO;

                /* std::vector<int> blockNums(blockNum, 0);
                  for (int i = start; i < end; ++i) {
                      blockNums[predict_block(kvs[i].first)]++;
                      //mean+=kvs[i].first/num_array_keys_;
                      //double pred =  seg_msg.slope* (kvs[i].first-minKey);
                      //double diff = pred - i+seg_msg.offset;
                      //loss += diff * diff;
                  }
                  double var = variance(blockNums, blockNum);
                  my_alloc::statistic.vars.push_back(var);
                  if (var > 38000)
                      int y = 0;*/
/*                my_alloc::allocMem(metas, blockNum);
                for (int i = 0; i < blockNum; ++i) {
                    metas[i] = new std::vector<LayerMeta>();
                }*/
                for (int i = start; i < end; i++) {
                    int idx = predict_block(kvs[i].first, blockNum);
                    auto block = blocks + idx;
                    auto blockHead = &block->head;
                    auto meta = metas + idx;
                    uint8_t fp = hashcode1B(kvs[i].first);
                    RETRY:
                    if (meta->num < kvsPerBlock) {
                        block->kv[meta->num] = kvs[i];
                        meta->fps[meta->num] = fp;
                        meta->num++;
                        // blockHead->num++;
                    } else {
                        ////my_alloc::statistic.maxoverflow++;
                        migrate(idx, metas, _slot);
                        //blockHead->num=0;
                        goto RETRY;
                    }
                }

            }
            // my_alloc::statistic.vars.push_back(slot->levelnum*1.0/blockNum);

            slot->setChildNodePtr(reinterpret_cast<uintptr_t>(this));
            //slot->setBufferPtr(reinterpret_cast<uint64_t>(metas));
            slot->size = blockNum;
            slot->minkey = minKey;
            slot->slope = slope_;
            Alloc::Persist(&slot, sizeof(innerSlot<KeyType>));
            return;


            /*        else if(num_array_keys_ <0)
                    {
                        int end=num_array_keys_+start;
                        linearReg_w_expanding(kvs,start,b,slope_,num_array_keys_, false);
        //        linearReg_w_expanding(keys, a, b, num_nonempty, fanout, true);
                        int last_k_id = 0;
                        KeyType last_key = kvs[start].first;
                        int pos = -1;
                        int last_pos = LR_PRED(b, slope_, last_key, num_array_keys_);
                        KeyType final_key = kvs[end-1].first;
        //    int final_pos = LR_PRED(a, b, final_key, fanout);
                        if (b < 0 || last_pos == LR_PRED(b, slope_,  final_key, num_array_keys_)) {
                            linearReg_w_expanding(kvs,start,b,slope_,num_array_keys_, true);

                        }
                        slope_=slope_/INIT_RATIO;
                        blockNum=num_array_keys_/INIT_RATIO/kvsPerBlock+1;
                        //minKey = kvs[start].first;
                        b=b/INIT_RATIO+1/INIT_RATIO;

                        blocks=reinterpret_cast<Block<KeyType,ValueType>*>(Alloc::PMZAlloc(blockNum*sizeof(Block<KeyType,ValueType>)));

                        slot->size=blockNum;
                        slot->b=b;
                        slot->slope=slope_;

                        std::vector<int> blockNums(blockNum,0);

                        //double loss=0;
                        for(int i=start;i<end;++i){
                            blockNums[predict_block(kvs[i].first)]++;
                            //double pred =  seg_msg.slope* (kvs[i].first-minKey);
                            //double diff = pred - i+seg_msg.offset;
                            //loss += diff * diff;
                        }
                        double var=variance(blockNums,blockNum);
                        my_alloc::statistic.vars.push_back(var);
                        if(var>38000)
                            int y=0;
                        //my_alloc::statistic.loss+=loss;
                        for(int i=start;i<end;i++){
                            int idx= predict_block(kvs[i].first);
                            auto block=blocks+idx;
                            auto blockHead = &block->head;
                            RETRY:
                            if(blockHead->num<kvsPerBlock) {
                                block->kv[blockHead->num]=kvs[i];
                                blockHead->num++;
                            }else{
                                my_alloc::statistic.maxoverflow++;
                                migrate(idx);
                                goto RETRY;
                            }
                        }
                    }*/

        }


        //计算数组Vector<int> arr的方差
        double variance(const std::vector<int> &arr, int num) {
            double sum = 0.0;
            double sum2 = 0.0;
            for (int i = 0; i < num; i++) {
                if (arr[i] == 0)
                    my_alloc::statistic.block0++;
                else if (arr[i] <= 15 && arr[i] != 0)
                    my_alloc::statistic.block15++;
                else my_alloc::statistic.blockmore15++;
                sum += arr[i];
                sum2 += arr[i] * arr[i];
            }
            return sum2 / num - sum * sum / num / num;
        }

        void var_mean(const std::vector<KV> &arr, int start, int num, double &mean, double &var) {
            double sum = 0.0;
            double sum2 = 0.0;
            uint64_t len = arr[start + num - 1].first - arr[start].first;
            for (int i = start; i < start + num; i++) {
                sum += arr[i].first * 1.0 / num;
            }
            mean = sum;
            for (int i = start; i < start + num; i++) {
                sum2 += (arr[i].first * 1.0 - mean) * (arr[i].first * 1.0 - mean) / num;
            }
            var = sqrt(sum2);
        }

        /*  inline void
          AddKV(const SegmentMessage<KeyType> &seg_msg, const std::vector<KV> &kvs, SegmentInf<KeyType> *&segInf) {
              //
              int kvsPerBlock = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
              num_array_keys_ = seg_msg.size;
              size_t end = num_array_keys_ + seg_msg.offset;
              blockNum = num_array_keys_ / INIT_RATIO / kvsPerBlock + 2;
              maxKey = seg_msg.key;
              minKey = kvs[seg_msg.offset].first;
              slope_ = seg_msg.slope / INIT_RATIO;
              segInf = new SegmentInf<KeyType>(maxKey, minKey, blockNum, 0, slope_, (uintptr_t) this);

              std::vector<int> blockNums(blockNum, 0);
              double loss = 0;
              for (int i = seg_msg.offset; i < end; ++i) {
                  blockNums[predict_block(kvs[i].first, blockNum - 1)]++;
                  double pred = seg_msg.slope * (kvs[i].first - minKey);
                  double diff = pred - i + seg_msg.offset;
                  loss += diff * diff;
              }
              double var = variance(blockNums);
              my_alloc::statistic.vars.push_back(var);
              my_alloc::statistic.loss += loss;

              blocks = reinterpret_cast<Block<KeyType, ValueType> *>(Alloc::PMZAlloc(
                      blockNum * sizeof(Block<KeyType, ValueType>)));
              for (int i = seg_msg.offset; i < end; i++) {
                  int idx = predict_block(kvs[i].first, blockNum - 1);
                  auto block = blocks + idx;
                  auto blockHead = &block->head;
                  RETRY:
                  if (blockHead->num < kvsPerBlock) {
                      block->kv[blockHead->num] = kvs[i];
                      blockHead->num++;
                  } else {
                      ////my_alloc::statistic.maxoverflow++;
                      migrate(idx);
                      goto RETRY;
                  }
              }
              Alloc::MemPoolAllocator::Persist(blocks, blockNum * sizeof(Block<KeyType, ValueType>));
              Alloc::MemPoolAllocator::Persist(this, sizeof(this));
          }

  */


        inline KV getKVAt(int idx, int pos) {
            return blocks[idx].kv[pos];
        }


        inline KeyType getKeyAt(int idx, int pos) {
            return blocks[idx].kv[pos].first;
        }

        inline ValueType getValueAt(int idx, int pos) {
            return blocks[idx].kv[pos].second;
        }

        inline void addBufferAt(int idx, uint64_t ptr, uint64_t fpPtr, blockMeta *bm) {
            blocks[idx].head.setBufferPtr(ptr);
            //blocks[idx].head.setFpPtr(fpPtr);
            bm[idx].setLayerPtr(ptr);
            bm[idx].setFPtr(fpPtr);
            Alloc::Persist(&blocks[idx].head, 8);
        }


        inline int insert(KeyType key, ValueType value, size_t idx, size_t pos) {
            blocks[idx].kv[pos].second = value;
            blocks[idx].kv[pos].first = key;
            Alloc::Persist(&blocks[idx].kv[pos], sizeof(KV));
            return 2;
        }

        inline void insertAt(KeyType key, ValueType value, size_t idx, size_t pos) {
            blocks[idx].kv[pos].second = value;
            blocks[idx].kv[pos].first = key;
            Alloc::Persist(&blocks[idx].kv[pos], sizeof(KV));
        }

        inline KV *getKVAddrAt(int idx, int pos) {
            return &blocks[idx].kv[pos];
        }

        inline void setEpochAt(int idx, uint16_t epoch) {
            blocks[idx].head.epoch = epoch;
            Alloc::Persist(&blocks[idx].head.epoch, sizeof(uint16_t));
        }

        inline void delKeyAt(int idx, int pos) {
            blocks[idx].kv[pos].first = FREE;
            Alloc::Persist(&blocks[idx].kv[pos].first, sizeof(KeyType));
        }

        inline void setValueAt(int idx, int pos, ValueType value) {
            blocks[idx].kv[pos].second = value;
            Alloc::Persist(&blocks[idx].kv[pos].second, sizeof(ValueType));
        }

        inline void deleteKVAt(int idx, int pos) {
            neiborInsert(FREE, reinterpret_cast<ValueType>(DELETED), idx, pos);
        }

        inline void setBucketAt(int idx, void *ptr) {
            blocks[idx].head.bucket = ptr;
            Alloc::Persist(&blocks[idx].head.bucket, sizeof(void *));
        }

        void neiborInsert(KeyType key, ValueType value, size_t idx, size_t pos) {
            // RETRY:
            blocks[idx].head.member.setBit(pos);
            if ((char *) &blocks[idx].head + 64 > (char *) &blocks[idx].kv[pos]) {
                blocks[idx].kv[pos].second = value;
                blocks[idx].kv[pos].first = key;
                Alloc::Persist(&blocks[idx].kv[pos], sizeof(KV));
            } else {
                Alloc::Persist(&blocks[idx].head, sizeof(BlockHeader));
            }
            blocks[idx].kv[pos].second = value;
            blocks[idx].kv[pos].first = key;
            Alloc::Persist(&blocks[idx].kv[pos], sizeof(KV));
        }

        inline void delKeyUnsafeAt(int idx, int pos) {
            blocks[idx].kv[pos].first = FREE;
        }

        inline void set_slope(float slope) {
            slope_ = slope;
            Alloc::Persist(&slope_, sizeof(float));
        }

        inline float get_slope() {
            return slope_;
        }

        inline void set_pre_segment(Segment<KeyType, ValueType> *pre) {
            pre_ = pre;
            Alloc::Persist(&pre_, sizeof(Segment<KeyType, ValueType> *));
        }

        inline void set_next_segment(Segment<KeyType, ValueType> *next) {
            next_ = next;
            Alloc::Persist(&next_, sizeof(Segment<KeyType, ValueType> *));
        }

        inline void switch_next_segment(Segment<KeyType, ValueType> *next) {
            if (flag == 0) {
                next1_ = next;
                Alloc::Persist(&next1_, sizeof(Segment<KeyType, ValueType> *));
                flag = 1;
                Alloc::Persist(&flag, sizeof(uint8_t));
            } else {
                next_ = next;
                Alloc::Persist(&next_, sizeof(Segment<KeyType, ValueType> *));
                flag = 1;
                Alloc::Persist(&flag, sizeof(uint8_t));
            }


        }

        inline Segment<KeyType, ValueType> *get_next_segment() {
            if (flag == 0)
                return next_;
            else
                return next1_;
        }

        inline Segment<KeyType, ValueType> *pre_segment() {
            return pre_;
        }

        inline Segment<KeyType, ValueType> *next_segment() {
            return next_;
        }

        inline uint32_t getBlockNum() { return blockNum; }

        inline uint32_t GetTotalKvNum() {
            return num_array_keys_ + num_buffers_keys_;
        }

        inline int numArray() {
            return num_array_keys_;
        }

        inline bool IsRetain(size_t avg_num_seg_keys) {
            // lazy retrain
            // Retrain when the number of sorted keys in buffer reaches a certain threshold to reduce sorting overhead.
/*            if (GetTotalKvNum() > avg_num_seg_keys * alpha_) {
//                std::cout << alpha_ * avg_num_seg_keys << std::endl;
                alpha_ = alpha_  * 2;
                return true;
            }*/
            return false;
        }

        inline KeyType back() {
            // return kvs_[num_array_keys_ - 1].first;
        }

        inline void setSlot(innerSlot<KeyType> *slot_) {
            slot = slot_;
            //Alloc::Persist(&slot, sizeof(innerSlot<KeyType> *));
        }

        inline innerSlot<KeyType> *getSlot() {
            return slot;
            //Alloc::Persist(&slot, sizeof(innerSlot<KeyType> *));
        }

        inline uint32_t getNum() {
            return blockNum;
        }

        inline KeyType get_minKey() {
            return minKey;
        }

        inline void set_metas(blockMeta *m) {
            metas = m;
        }


    private:
        // funtype functions[2];
        KeyType maxKey;
        KeyType minKey;
        float slope_;
        uint32_t blockNum;
        blockMeta *metas;
        Segment<KeyType, ValueType> *pre_;
        Segment<KeyType, ValueType> *next_;
        Segment<KeyType, ValueType> *next1_;
        innerSlot<KeyType> *slot;
        uint32_t num_array_keys_;
        uint32_t num_buffers_keys_;
        uint32_t version;
        uint8_t flag;
        uint8_t unused[179];
        Block<KeyType, ValueType> blocks[0];
        //BufferNode<KeyType,ValueType,STASH_SIZE> * stash; //
    };
}

#endif //ART_TEST_SEGMENT_H
