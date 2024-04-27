
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
/**
 * 构造函数 - 创建一个指定层级的Layer对象。
 *
level 层级，表示当前Layer在结构中的位置。
fpPtr 指向FPEntrys的指针，通过这个指针可以访问到Layer中的FPEntrys数组。
lastFPEntrys 指向上一层级的FPEntrys数组的指针，如果当前层级是第0层，则为NULL。
_slot 指向内部槽的指针，用于维护当前Layer所在结构的元数据。
 */
        Layer(const size_t &level, uint64_t &fpPtr, EntryFp *lastFPEntrys, innerSlot<KeyType> *_slot) : level(level),
                // 根据层级计算当前Layer应包含的entry数量，并分配相应内存
            size_t size = 1 << (INIT_ENTRY_BIT + FANOUT * level);
            entrys = reinterpret_cast<EntryType *>(Alloc::PMZAlloc(size * sizeof(EntryType)));
            FPEntrys = NULL;
            my_alloc::allocMem(FPEntrys, (int) size);
            fpPtr = reinterpret_cast<uint64_t>(FPEntrys);
            for (int i = 0; i < size; ++i) {
                FPEntrys[i].entryPtr = reinterpret_cast<uintptr_t>(&entrys[i]);
            }
        // 如果当前层级不是第一层，则将上一层的每个entry的nextLevel指向当前的FPEntrys
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
       //在一个32位无符号整数的指定范围内反转位。
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
/**
 * 遍历所有条目获取所有的键值对并将其存储到提供的map中
 *
 deleted 用于统计被标记为删除的键值对数量。
 result  用于存储检索到的键值对。
 */
        void getAllKVs(int &deleted, std::map<KeyType, ValueType> &result) {
    // 根据当前等级计算总条目数
            int entryNum = 1 << (INIT_ENTRY_BIT + FANOUT * level);
            for (int i = entryNum - 1; i >= 0; --i) {
                // 跳过空条目
                if (metas[i].size == 0)
                    continue;
                // 计算有效桶的数量
                int n = metas[i].size % BUCKET_LEN == 0 ? metas[i].size / BUCKET_LEN - 1 : metas[i].size / BUCKET_LEN;
                // 遍历每个桶
                for (int j = n; j >= 0; --j) {
                    int k = (j < n ? BUCKET_LEN - 1 : metas[i].size % BUCKET_LEN - 1);
                    for (; k >= 0; --k) {
                        KV *kv_ = getKVAt(i, j, k);
                        // 忽略空键值对
                        if (entrys[i].buckets[j].kvs[k].first != FREE) {
                            // 尝试将键值对插入结果map
                            result[entrys[i].buckets[j].kvs[k].first] = entrys[i].buckets[j].kvs[k].second;
                            const auto [ret, success] = result.insert(
                                    std::make_pair<KeyType, ValueType>(std::move(kv_->first), std::move(kv_->second)));
                            // 如果插入成功且值被标记为删除，则增加删除计数
                            if (success && kv_->second == DELETED) {
                                ++deleted;
                            }
                        }
                    }
                }
            }
        }
/**
 * 在溢出缓存中进行范围扫描
 *
key 扫描的起始键值
num期望返回的结果数量
deleted扫描过程中被标记为删除的条目数量
results存储扫描结果的map
 */
        void rangeScan(KeyType key, int num, int &deleted, std::map<KeyType, ValueType> &results) {
    // 初始化结果集大小和上限边界
            int result_size = results.size();
            uint64_t upper_bound = UINT64_MAX;
    // 如果结果集不为空，找到最近的非删除条目的上界
            if (result_size > 0) {
                auto it = prev(results.end());
                while (it != results.begin() && it->second == DELETED) {
                    --it;
                }
                if (it->second != DELETED)
                    upper_bound = it->first;
            }
    // 计算桶的个数并遍历每个桶
            int entryNum = 1 << (INIT_ENTRY_BIT + FANOUT * level);
            for (int i = entryNum - 1; i >= 0; --i) {
                if (metas[i].size == 0)
                    continue;
                int n = metas[i].size % BUCKET_LEN == 0 ? metas[i].size / BUCKET_LEN - 1 : metas[i].size / BUCKET_LEN;
                for (int j = n; j >= 0; --j) {
                    int k = (j < n ? BUCKET_LEN - 1 : metas[i].size % BUCKET_LEN - 1);
                    for (; k >= 0; --k) {
                        KV *kv_ = getKVAt(i, j, k);
                        // 如果条目未被标记为FREE，且键值在范围内，则尝试插入结果集
                        if (kv_->first != FREE && kv_->first >= key &&
                            (result_size - deleted < num || kv_->first < upper_bound)) {
                            const auto [ret, success] = results.insert(
                                    std::make_pair<KeyType, ValueType>(std::move(kv_->first), std::move(kv_->second)));
                            if (success) {
                                ++result_size;
                                if (kv_->second == DELETED)
                                    ++deleted;
                            }
                            // 如果结果集超过期望数量，删除多余的元素
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
/**
 * 将指定条目的数据批量合并到下一层级中。
 */
        void mergeToNextLevel(int entryIdx, innerSlot<KeyType> *_slot) {
    // 如果下一层级尚未初始化，则创建一个新的目录层级
            if (nextlevel == nullptr) {
                uint64_t fpPtr = 0;
                nextlevel = new DirType(level + 1, fpPtr, FPEntrys, _slot);
                Alloc::Persist(&nextlevel, sizeof(DirType *));
            }
    // 用于存储重新哈希后的键值对数组
            KV *kvs[BUCKET_PER_ENTRY * (BUCKET_LEN) * (1 << FANOUT)] = {nullptr};
            uint32_t rehashed = 0;
            uint32_t bucketIdx = 0;
            uint32_t sizes[1 << FANOUT] = {0};
    // 遍历当前条目下的所有桶，进行重新哈希
            while (rehashed < metas[entryIdx].size && bucketIdx < BUCKET_PER_ENTRY) {
                int toRehash = std::min((uint32_t) BUCKET_LEN, metas[entryIdx].size - rehashed);
                rehash(nextlevel->level, entryIdx, bucketIdx, toRehash, kvs, sizes);
                rehashed += toRehash;
                ++bucketIdx;
            }
    // 将重新哈希后的键值对批量插入到下一层级
            bulkInsertToNextLevel(kvs, sizes, level + 1, _slot);
            FPEntrys[entryIdx].clear();
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
/**
 * 将一批键值对批量插入到下一层级中。
 kvs 指向键值对数组的指针数组。
 sizes 键值对数组对应每个子数组的大小。
 targetlevel 目标层级，即要插入的下一层级。
 _slot 内部节点插槽位指针。
 */
        void bulkInsertToNextLevel(KV **kvs, uint32_t *sizes,
                                   size_t targetlevel, innerSlot<KeyType> *_slot) {
    // 计算每个子数组中键值对的理论数量
            uint32_t kvnum = (BUCKET_LEN) * BUCKET_PER_ENTRY;
    // 计算当前层级的每个条目的大小
            uint32_t entrySize = 1 << (FANOUT);
            for (uint32_t i = 0; i < entrySize; ++i) {
                if (sizes[i] == 0) {
                    continue;
                }
                // 计算当前子数组的起始索引
                uint32_t start = i * kvnum;
                // 获取当前子数组对应的键值对
                auto kv = kvs[start];
                // 计算在下一层级中该键值对应该插入的位置
                uint32_t entryIdx = getEntryIdx(nextlevel->level, kv->first);
                // 获取对应位置的元数据
                EMetaType *meta = nextlevel->metas + entryIdx;
                uint32_t size_ = meta->size;
                // 如果合并后的大小超过限制，则触发在下一层级的合并操作
                if (size_ + sizes[i] > BUCKET_PER_ENTRY * BUCKET_LEN) {
                    nextlevel->mergeToNextLevel(entryIdx, _slot);
                }
                // 执行批量插入操作，并返回插入的键值对数量
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
/**
 * 批量插入键值对到下一个级别中的指定条目。
 */
        int
        bulkInsert(DirType *dir, uint32_t nextLevelEntryIdx, uint32_t size, uint32_t start, KV **kvs,
                   size_t targetlevel) {
            EntryType *entry = dir->entrys + nextLevelEntryIdx; // 获取指定条目的指针
            EMetaType *meta = dir->metas + nextLevelEntryIdx;
            int n = getFreeBucketIdx(meta->size);
            if (n == -1) {
                // 所有桶均已满
                return 0;
            }
            uint32_t elemsInserted = 0; // 已插入元素计数
            while (n < BUCKET_PER_ENTRY && elemsInserted < size) {
                BucketType *bucket = entry->buckets + n; // 获取当前操作的桶
                uint32_t bucketSize = getSizeOfBucket(meta->size, n); // 获取当前桶已使用的大小
                auto fp_bucket = dir->FPEntrys[nextLevelEntryIdx].fps + n * BUCKET_LEN;// 获取桶的指纹数组位置
                uint32_t elemsToInsert = std::min(BUCKET_LEN - bucketSize, size - elemsInserted); // 计算可插入元素数量
                if (elemsToInsert > 0) {
                    // 将键值对批量插入到桶中
                    for (uint32_t i = 0; i < elemsToInsert; ++i) {
                        uint32_t offset = start + elemsInserted + i;
                        _mm_stream_si64((long long *) (&bucket->kvs[bucketSize + i].first),
                                        (long long) kvs[offset]->first);
                        _mm_stream_si64((long long *) (&bucket->kvs[bucketSize + i].second),
                                        (long long) kvs[offset]->second);
                    }
                    // 为插入的键值对计算并设置指纹
                    for (uint32_t i = 0; i < elemsToInsert; ++i) {
                        uint32_t offset = start + elemsInserted + i;
                        unsigned char fp = hashcode1B<KeyType>(kvs[offset]->first);
                        fp_bucket[bucketSize + i] = fp;
                    }
                    elemsInserted += elemsToInsert; // 更新已插入元素计数
                }
                ++n; // 移动到下一个桶
            }
            meta->size += elemsInserted; // 更新元数据中entry中有效kv数
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
//重建块元数据
  void rebuildBlockMeta(int idx1, uint8_t version, blockMeta *bm1) {
            // 计算每个块中可以存储的KV对数量
            uint8_t kvPerBlock = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
            uint8_t num = 0;
            auto bm=metas;
            // 遍历所有块
            for(int idx=0;idx<blockNum;idx++){
                auto ptr = blocks[idx].head.ptr;
                bm[idx].setLayerPtr(ptr);
                bm[idx].v = version;
                // 遍历块中的KV对，统计有效KV对数量并更新元数据
                for (uint8_t i = 0; i < kvPerBlock; ++i) {
                    if (blocks[idx].kv[i].first != FREE) {
                        uint8_t fp = hashcode1B<KeyType>(blocks[idx].kv[i].first);
                        bm[idx].fps[i] = fp;
                        bm[idx].num++;
                    } else {
                        break;
                    }
                }
                EntryFp *lastFPEntrys = NULL;
                auto layer = reinterpret_cast< DirType *> (ptr);
                int level = 0;
                // 遍历每个层级的entry数组，构建缓存元数据
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
                        bm[idx].setFPtr(reinterpret_cast<uint64_t>(layer->FPEntrys));
                    } else if (level > 0 && lastFPEntrys != NULL) {
                        int lastSize = 1 << (INIT_ENTRY_BIT + FANOUT * (level - 1));
                        for (int i = 0; i < lastSize; ++i) {
                            lastFPEntrys[i].nextLevel = layer->FPEntrys;
                        }

                    }
                    lastFPEntrys = layer->FPEntrys;
                    layer = layer->nextlevel;
                    level++;
                }
            }

        }

        inline bool key_equal(const KeyType &a, const KeyType &b) const {
            return !(a < b) && !(b < a);
        }

        /**
  * 函数名称: upsert
  * 功能描述: 在指定的槽(_slot)内插入或更新键值对(key, value)。如果键已存在，则更新其值；如果不存在，则插入新的键值对。
  * 参数列表:
  *    key - 要插入或更新的键；
  *    value - 要插入或更新的值；
  *    _slot - 指向内部槽的指针，用于确定键值对存储的位置。
  * 返回值:
  *    OpStatus - 操作状态枚举值，表示操作的成功与否或具体失败原因。
  */
        OpStatus
        upsert(const KeyType &key, const ValueType &value,  innerSlot<KeyType> *_slot) {
            // 计算每个块能容纳的键值对数量
            uint8_t kvPerBlock = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
            // 获取块元数据指针
            auto bm=reinterpret_cast<blockMeta * >(_slot->getBufferPtr());
            // 根据键计算出块的位置
            auto block=(key-_slot->minkey) * _slot->slope /kvPerBlock;
            // 确保索引在有效范围内，防止溢出
            int idx=block<0?0:(block>=_slot->size?_slot->size-1:block);
            // 检查索引是否有效和槽是否为空，不合法则返回失败
            if (idx == -1 || _slot == NULL)
                return OpStatus::Failed;
            // 检查读锁是否可用，不可用则返回锁定状态
            if (!_slot->check_read_lock())
                return OpStatus::lock;
            // 检查是否需要进行叶节点扩展
            if (_slot->levelnum / _slot->size > THRESHOLD)
                return OpStatus::SMO;
            // 获取块元数据的指针
            auto meta = bm + idx;
            // 检查写锁是否可用，不可用则返回锁定状态
            if (!_slot->check_write_lock())
                return OpStatus::lock;
            // 获取元数据中block的锁
            meta->get_lock();
            // 如果元数据版本不匹配，则重建块元数据
            if (meta->v != _slot->v) {
                rebuildBlockMeta(idx, _slot->v, meta);
            }
            RETRY:
            // 查找空闲的位置
            uint8_t pos = meta->num;
            // 如果没有空闲位置，则将block的数据批量写入溢出缓存并重试
            if (pos >= kvPerBlock) {
                migrate(idx, bm, _slot);
                goto RETRY;
            }
            // 计算键的指纹
            uint8_t fp = hashcode1B<KeyType>(key);
            // 遍历已有的键值对，查找是否存在相同的键
            for (uint8_t i = 0; i < pos; ++i) {
                if (meta->fps[i] == fp) {
                    if (key_equal(blocks[idx].kv[i].first, key)) {
                        // 如果找到相同的键，则更新其值并返回更新状态
                        blocks[idx].kv[i].second = value;
                        Alloc::Persist(&blocks[idx].kv[i].second, sizeof(ValueType));
                        meta->release_lock();
                        return OpStatus::Updated;
                    }
                }
            }
            // 如果没有找到相同的键，则插入新的键值对
            blocks[idx].kv[pos].second = value;
            blocks[idx].kv[pos].first = key;
            meta->fps[pos] = fp;
            Alloc::Persist(&blocks[idx].kv[pos], sizeof(KV));
            meta->addNum();
            meta->release_lock();
            // 插入成功，返回插入状态
            return OpStatus::Inserted;
        }
/**
 * 查找给定键对应的值。
 *
 *key 要查找的键。
 * value 如果找到键，则通过引用返回对应的值。
 *_slot 存储叶节点信息的内部槽位的指针。
 * @return 返回操作的状态，可以是成功找到、锁冲突、需要重试或未找到。
 */
        OpStatus find(const KeyType &key, ValueType &value, innerSlot<KeyType> *_slot) {
            // 获取块元数据指针
            auto bm=reinterpret_cast<blockMeta * >(_slot->getBufferPtr());
            // 计算块的位置
            auto block=(key-_slot->minkey) * _slot->slope /15;
            // 确保索引在有效范围内
            int idx=block<0?0:(block>=_slot->size?_slot->size-1:block);
            // 如果索引无效，则返回失败状态
            if (idx == -1)
                return OpStatus::Failed;
            // 检查读锁是否可用
            if (!_slot->check_read_lock()) {
                return OpStatus::lock;
            }
            // 获取块元数据，并检查版本是否一致
            auto meta = bm + idx;
            if (__glibc_unlikely(meta->v != _slot->v)) {
                // 尝试获取锁，如果失败则返回重试状态
                if (!meta->try_get_lock()) {
                    return OpStatus::Retry;
                }
                // 重建块元数据
                rebuildBlockMeta(idx, _slot->v, meta);
                // 释放锁而不改变状态
                meta->release_lock_without_change();
            }
            // 加载元数据的版本和关键字数量
            uint32_t old_version = __atomic_load_n(&meta->version_lock, __ATOMIC_ACQUIRE);
            uint8_t pos = __atomic_load_n(&meta->num, __ATOMIC_ACQUIRE);
            // 如果版本标记为锁定，则返回锁定状态
            if (old_version & lockSet) {
                return OpStatus::lock;
            }
            // 计算每个块内关键字的数量
            auto kvPerBlock = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
            // 计算键的指纹
            uint8_t fp = hashcode1B<KeyType>(key);
            // 在当前块内查找键值对
            for (int i = 0; i < kvPerBlock; ++i) {
                if (meta->fps[i] == fp && key_equal(blocks[idx].kv[i].first, key)) {
                    value = blocks[idx].kv[i].second;
                    // 检查block锁版本或pos是否发生变化，若已变化则返回重试状态
                    if (meta->check_version_change(old_version) || meta->check_size_change(pos))
                        return OpStatus::Retry;
                    return OpStatus::Find;
                }
            }
            // 在溢出缓存中查找键值对
            auto ret1 = false;
            uint64_t buffer = meta->getLayerPtr();
            if (__glibc_likely(buffer!=0)) {
                ret1 = searchInDir(key, value, reinterpret_cast<EntryFp *>(meta->getFPtr()));
                // 检查版本或pos是否发生变化，若已变化则返回重试状态
                if (meta->check_version_change(old_version) || meta->check_size_change(pos))
                    return OpStatus::Retry;
            }
            // 根据查找结果返回相应的状态
            return ret1 ? OpStatus::Find : OpStatus::NotFound;
        }

//键值对的删除用插入DELETE标记的方式进行
        inline OpStatus erase(const KeyType &key,  innerSlot<KeyType> *_slot) {
            return upsert(key, reinterpret_cast<ValueType>(DELETED),  _slot);
        }
//遍历block
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
            int bucketIdx = n / BUCKET_LEN;
            BucketType *bucket = &_entry->buckets[bucketIdx];
            if (bucket->kvs[n % BUCKET_LEN].first == key) {
                value = bucket->kvs[n % BUCKET_LEN].second;
                return true;
            }
            return false;
        }

/**
 * 在溢出缓存中逐层查找关键字
 *
 * key 要查找的关键字
 * value 用于存储找到的关键字对应的价值，如果找到的是已删除的条目，则不修改value的值
 *fp_entrys 缓存元数据指针
 * @return 如果找到且未被删除，则返回true；如果找到但已删除，则返回false；如果未找到，则返回false。
 */
        bool searchInDir(KeyType key, ValueType &value, EntryFp *fp_entrys) {
            // 使用hashcode1B函数计算关键字的哈希值
            unsigned char fp = hashcode1B<KeyType>(key);
            __m128i target = _mm_set1_epi8(fp);
            int level = 0; // 当前所在的层级
            // 当前目录条目非空时继续查找
            while (fp_entrys != NULL) {
                // 根据层级和关键字计算目录条目的索引
                int entryIdex = DirType::getEntryIdx(level, key);
                // 遍历每个桶中的条目
                for (int i = BUCKET_PER_ENTRY - 1; i >= 0; i--) {
                    __m128i data = _mm_loadu_si128((__m128i *) &fp_entrys[entryIdex].fps[i * 16]); // 加载数据
                    __m128i compare = _mm_cmpeq_epi8(data, target); // 比较数据
                    int mask = _mm_movemask_epi8(compare); // 获取比较结果的掩码

                    // 当比较结果中存在匹配时
                    while (mask != 0) {
                        int pos1 = __builtin_ctz(mask); // 找到掩码中最低位的1的位置
                        int position = i * 16 + pos1; // 计算匹配位置
                        // 在当前层中查找具体的条目
                        auto ret = findInLayer(key, value, position,
                                               reinterpret_cast<EntryType *>(fp_entrys[entryIdex].entryPtr));
                        if (ret) {
                            // 如果找到的条目未被删除，则返回true
                            if (value != DELETED)
                                return true;
                            else
                                return false; // 如果找到的条目已被删除，则返回false
                        }
                        // 更新掩码，移除已经检查过的位
                        mask = mask & ~(1 << pos1);
                    }
                }
                // 移动到下一层级继续查找
                level++;
                fp_entrys = fp_entrys[entryIdex].nextLevel;
            }
            // 如果未在任何层级找到关键字，则返回false
            return false;
        }
        /**
        * 将指定索引的block中的数据批量写入溢出缓存。
        * idx 指定的block索引。
        * bm 指向blockMeta的指针，用于管理block的元数据。
        *_slot 指向innerSlot的指针，其中存储了叶节点信息。
        */
        void migrate(const int &idx, blockMeta *bm, innerSlot<KeyType> *_slot) {
            // 获取指定索引block的头部信息
            auto blockHead = &blocks[idx].head;
            // 获取溢出缓存指针
            auto dir = reinterpret_cast<DirType *> (blockHead->getBufferPtr());
            // 如果dir==null创建新的层
            newDir(dir, idx, 0, bm, _slot);

            // 初始化用于存储KV对指针的数组和记录每个entry中KV数量的数组
            KV *kvs[1 * BUCKET_LEN * (1 << INIT_ENTRY_BIT)] = {nullptr};
            uint32_t sizes[1 << (INIT_ENTRY_BIT)] = {0};

            // 计算当前block中KV对的数量
            uint32_t num = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
            // 计算上一层entry中的KV对数量
            uint32_t kvnum = 1 * BLOCK_SIZE / sizeof(KV);

            // 遍历当前block中的所有KV对，将其按照目录索引存储到kvs数组中，并记录每个entry的KV数量
            for (int i = 0; i < num; ++i) {
                uint32_t newEntryIdx = DirType::getEntryIdx(0, getKeyAt(idx, i));
                kvs[newEntryIdx * kvnum + sizes[newEntryIdx]] = blocks[idx].kv + i;
                sizes[newEntryIdx]++;
            }

            // 处理每个entry，将KV对批量插入到entry中
            int entrySize = 1 << INIT_ENTRY_BIT;
            for (int i = 0; i < entrySize; ++i) {
                if (sizes[i] == 0) {
                    continue;
                }
                uint32_t start = i * kvnum;
                auto key = kvs[start]->first;
                uint32_t entryIdx = DirType::getEntryIdx(0, key);
                EntryMeta<KeyType> *meta = dir->metas + entryIdx;

                // 如果插入后的大小超过单个entry的容量，则将数据合并到下一级
                if (meta->size + sizes[i] > BUCKET_PER_ENTRY * BUCKET_LEN) {
                    dir->mergeToNextLevel(entryIdx, _slot);
                }

                // 批量插入KV对到entry中
                int retNum = dir->bulkInsert(dir, entryIdx, sizes[i], start, kvs, 0);
                // 确保插入的KV对数量与预期一致
                assert (retNum == sizes[i]);
            }

            // 从block中删除已插入的KV对
            for (int i = 0; i < num; ++i) {
                delKeyUnsafeAt(idx, i);
            }

            // 持久化block数据
            Alloc::Persist(&blocks[idx], BLOCK_SIZE);
            // 清理block元数据
            bm[idx].clearNum();
            bm[idx].clear();
        }

/**
 * 更新结果集合。
 * 该函数遍历 `kvs` 中的键值对，并根据 `lower_bound`、`num` 和 `deleted` 的约束，
 * 将符合条件的键值对插入到 `results` 中。当 `results` 的大小超过 `num` 时，
 * 删除最老的（即键值对比较大的）元素，以保持 `results` 的大小不超过 `num`。
 *
 * lower_bound 插入键值对的下界。
 *  kvs 键值对集合，将从中选择并插入到 `results` 中。
 * num `results` 集合允许的最大大小。
 * deleted 已经从 `results` 中删除的元素数量。
 * results 用于存储结果的键值对集合。
 */
        void updateResults(const KeyType &lower_bound, std::vector<KV> &kvs, const int &num, const int &deleted,
                           std::map<KeyType, ValueType> &results) {
            int result_size = results.size(); // 当前结果集合的大小
            uint64_t upper_bound = UINT64_MAX; // 插入键值对的上界，默认为最大无符号64位整数

            // 如果结果集合不为空，查找最后一个非删除状态的元素，更新上界
            if (result_size > 0) {
                auto it = prev(results.end());
                while (it != results.begin() && it->second == DELETED) {
                    --it;
                }
                if (it->second != DELETED)
                    upper_bound = it->first;
            }

            // 遍历 kvs，将符合条件的键值对插入 results
            for (int i = 0; i < kvs.size(); ++i) {
                // 判断是否满足插入条件
                if (kvs[i].first >= lower_bound &&
                    (result_size - deleted < num || kvs[i].first < upper_bound)) {
                    const auto [it, success] = results.insert(
                            std::make_pair<KeyType, ValueType>(std::move(kvs[i].first), std::move(kvs[i].second)));
                    if (success) { // 成功插入新元素
                        ++result_size;
                    }
                    // 如果结果集合超过指定大小，删除最老的元素
                    if (result_size - deleted > num) {
                        results.erase(prev(results.end()));
                        upper_bound = prev(results.end())->first;
                        --result_size;
                    }
                }
            }
        }

/**
 * 从blocks和对应溢出缓存中结构中获取所有键值对，并将它们存储到一个给定的map中。
 *
 * results 引用传递，用于存储检索到的键值对的map。
 *blocknum 需要遍历的块的数量。
 * return 返回成功添加到results中的键值对的数量。
 */
        int getAllKVs(std::map<KeyType, ValueType> &results, const uint32_t &blocknum) {
            // 计算每个块中键值对的最大数量
            int blockSize = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
            // 用于追踪被标记为删除的键值对的数量
            int deleted = 0;

            // 遍历所有指定的块
            for (int idx = 0; idx < blocknum; ++idx) {
                // 遍历当前块中的所有键值对
                for (int i = 0; i < blockSize; ++i) {
                    // 获取当前位置的键值对
                    auto kv = getKVAt(idx, i);
                    // 如果键值对未被标记为FREE，则尝试将其添加到results中
                    if (kv.first != FREE) {
                        const auto [it, success] = results.insert(
                                std::make_pair<KeyType, ValueType>(std::move(kv.first), std::move(kv.second)));
                        // 如果插入成功且值被标记为DELETED，则增加deleted计数
                        (success && kv.second == DELETED) ? deleted++ : deleted;
                    }
                }
                // 检查是否存在溢出缓存
                if (blocks[idx].head.ptr != 0) {
                    auto ptr = reinterpret_cast<DirType *>(blocks[idx].head.ptr);
                    while (ptr != nullptr) {
                        // 递归调用，从下层块中获取键值对
                        ptr->getAllKVs(deleted, results);
                        ptr = ptr->nextlevel; // 移动到下一层
                    }
                }
            }
            // 返回最终结果，即results中的键值对数量
            return results.size();

        }

/**
 * 对给定的键值范围进行扫描，并将结果存储在一个map中。
 *
 *key 用于范围扫描的键值起点。
 * num期望扫描返回的元素数量。
 * deleted 已删除元素的计数器。
 *results 存储扫描结果的map。
 *idx 当前正在扫描的块的位置。
 * blocknum 需要扫描的块的总数。
 * return 实际扫描返回的元素数量。
 */
        int
        rangeScan(const KeyType &key, const int &num, int &deleted, std::map<KeyType, ValueType> &results, int idx,
                  const int &blocknum) {
            // 如果索引为-1，表示没有更多的块需要扫描，直接返回0
            if (idx == -1)
                return 0;

            // 在当前索引位置开始扫描，并更新相关参数
            int ret = scanAtBlock(idx, key, num, deleted, results, blocknum);
            idx++;

            // 继续扫描直到达到指定的元素数量或扫描完所有块
            while (ret < num&&idx < blocknum  ) {
                ret = scanAtBlock(idx, key, num, deleted, results, blocknum);
                idx++;
            }

            // 返回实际扫描到的元素数量
            return ret;
        }

/**
 * 在指定block中扫描匹配给定键值的记录。
 *
 *  idx 指定块的索引。
 * key 要匹配的键值。
 * num 需要返回的匹配项数量。
 * deleted 用于记录已删除的项的数量。
 *results 存储扫描结果的map。
 * blocknum 要扫描的块的数量。
 * return 返回实际返回的未删除项的数量。
 */
        int
        scanAtBlock(const int &idx, const KeyType &key, const int &num, int &deleted,
                    std::map<KeyType, ValueType> &results, const int &blocknum) {
            std::vector<KV> kvs; // 用于临时存储扫描到的项。
            int blockSize = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV); // 计算每个块中实际存储KV对的数量。

            // 从块中扫描匹配的KV对。
            if (blocknum > 0) {
                for (int i = 0; i < blockSize; ++i) {
                    auto kv = getKVAt(idx, i); // 获取当前位置的KV对。
                    if (kv.first != FREE && kv.first >= key) { // 筛选符合要求的KV对。
                        kvs.push_back(kv); // 加入到临时存储中。
                        kv.second == DELETED ? deleted++ : deleted; // 如果是已删除项，则计数。
                    }
                }
            }

            // 更新结果集。
            updateResults(key, kvs, num, deleted, results);
            // 遍历并扫描下一级块。
            if (blocks[idx].head.ptr != 0) {
                auto ptr = reinterpret_cast<DirType *>(blocks[idx].head.ptr);
                while (ptr != nullptr) {
                    ptr->rangeScan(key, num, deleted, results); // 对下一级块进行范围扫描。
                    ptr = ptr->nextlevel; // 移动到下一个块。
                }
            }
            return results.size() - deleted; // 返回未删除的项的数量。
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
/**
 * 对叶节点进行扩展操作
 */
        void resize(const std::map<KeyType, ValueType> kvs, innerSlot<KeyType> *_slot, Segment<KeyType, ValueType> *pre,
                    Segment<KeyType, ValueType> *next) {
            // 计算每个块能容纳的键值对数量
            int kvsPerBlock = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
            pre_ = pre;
            next_ = next;
            auto begin = kvs.begin();
            while (begin != kvs.end() && begin->second == DELETED) {
                begin++;
            }
            assert(begin != kvs.end());
            auto end = prev(kvs.end());
            // 找到最后一个非删除的键值对
            while (end != kvs.begin() && end->second == DELETED) {
                --end;
            }
            assert(end != kvs.begin());
            // 重新计算模型斜率
            slope_ = num_array_keys_ * 1.0 / (end->first - begin->first);
            minKey = begin->first;
            slot = _slot;
            my_alloc::allocMem(metas, blockNum);
            for (int i = 0; i < blockNum; ++i) {
                metas[i].v = _slot->v;
                blocks[i].head.v = _slot->v;
            }
            slot->setBufferPtr(reinterpret_cast<uint64_t>(metas));
            // 遍历键值对并分配到相应的块中
            while (begin != kvs.end()) {
                if(begin->second==DELETED){
                    begin++;
                    continue;
                }
                int idx = predict_block(begin->first, blockNum);
                auto block = blocks + idx;
                auto blockHead = &block->head;
                auto meta = metas + idx;
                uint8_t fp = hashcode1B(begin->first);
                RETRY:
                // 如果块还有空间，则将键值对添加到块中
                if (meta->num < kvsPerBlock) {
                    block->kv[meta->num].first = begin->first;
                    block->kv[meta->num].second = begin->second;
                    meta->fps[meta->num] = fp;
                    meta->num++;
                } else
                    //批量写入下一层
                    migrate(idx, metas, _slot);
                    goto RETRY;
                }
                begin++;
            }

            return;
        }
//执行批量加载操作
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
                int end = num_array_keys_ + start;
                slope_ = slope_ / INIT_RATIO;
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
                    } else {
                        migrate(idx, metas, _slot);
                        goto RETRY;
                    }
                }

            }
            slot->setChildNodePtr(reinterpret_cast<uintptr_t>(this));
            slot->size = blockNum;
            slot->minkey = minKey;
            slot->slope = slope_;
            Alloc::Persist(&slot, sizeof(innerSlot<KeyType>));
            return;
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
//根据当前的标志位（flag）进行兄弟叶节点指针切换
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
