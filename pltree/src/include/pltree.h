//
// Created by mcooiedo on 22-12-4.
//

#ifndef PLTREE_H
#define PLTREE_H

#include <algorithm>
#include <cmath>
#include <iostream>

#include "segment.h"
#include "allocator.h"
//#include <mcheck.h>
//#define BTREE
#define  ART
#define ROOT_SIZE 7
#define  LEVEL_SIZE 20

namespace plt {



    template<typename KeyType, typename ValueType, typename Alloc=my_alloc::MemPoolAllocator>//PMAllocator>
    class PLTree {
        typedef struct KV<KeyType, ValueType> KV;
        //typedef  BufferNode<KeyType,ValueType> BufferNodeType;
        typedef Segment<KeyType, ValueType> SegType;
        typedef std::vector<std::pair<KeyType *, int>> TreeLayoutType;

    public:
        PLTree(bool checkrecover=true)
               {
                   // 如果是恢复状态且checkrecover为true，则进行版本递增和恢复检查
            if (Alloc::isRecover()&&checkrecover) {
                version=version+1==0?1:version+1;
                recoveryCheck();
            } else {
                // 默认初始化设置
                version=1;
                min_key_=std::numeric_limits<KeyType>::max();
                max_key_=std::numeric_limits<KeyType>::min();
                max_error_=256;
                segments_head_=nullptr;
                segments_tail_=nullptr;
            
            }
        }


        PLTree(const PLTree &) = delete;

        PLTree &operator=(const PLTree &) = delete;

        ~PLTree() {
            // todo free segments
            // ...
            auto cur_seg = segments_head_;
            while (cur_seg) {
                auto next_seg = cur_seg->get_next_segment();
                delete cur_seg;
                cur_seg = next_seg;
            }
            SegType::segment_allocated_byte = 0;
        }

/**
 * 使用FSW算法自底向上对数据进行分段，构建索引布局。
  error 允许的索引构建误差，默认值为256。
 */
        void setTreeLayout(const std::vector<KV> &kvs,size_t error=256) {
            // 确定键的最小值和最大值
            KeyType min_key = kvs[0].first;
            KeyType max_key = kvs.back().first;

            // 使用FSW算法的构建器初始化
            plt::Builder<KeyType> asb(min_key, max_key, error/*,alloc*/);

            // 分配内存以存储键
            KeyType *keys = NULL;
            my_alloc::align_zalloc((void **) &keys, kvs.size() * sizeof(KeyType));

            // 向构建器添加键，并存储到keys数组中
            int keyNum = 0;
            for (const auto &kv: kvs) {
                asb.AddKey(kv.first);
                keys[keyNum++] = kv.first;
            }
            asb.Finalize(); // 完成第一阶段的布局构建

            // 获取第一阶段构建的段信息
            auto &seg_message = asb.get_segments_message();

            // 将第一阶段的键信息存储到treeLayout中
            treeLayout<KeyType>.push_back(std::make_pair(keys, keyNum));

            // 循环，直到满足根节点大小要求
            while (seg_message.size() + 1 > ROOT_SIZE) {
                // 为每一级布局分配键存储空间
                KeyType *keys = NULL;
                my_alloc::align_zalloc((void **) &keys, (seg_message.size() + 1) * sizeof(KeyType));
                int keyNum = 0;

                // 使用段信息构建下一级的索引布局
                plt::Builder<KeyType> asb(min_key, seg_message.back().key, error);
                asb.AddKey(min_key);
                keys[keyNum++] = min_key;
                for (const auto &item: seg_message) {
                    asb.AddKey(item.key);
                    keys[keyNum++] = item.key;
                }
                asb.Finalize(); // 完成下一级的布局构建

                // 更新段信息以进行下一轮构建
                seg_message = asb.get_segments_message();

                // 将构建好的键信息存储到treeLayout中
                treeLayout<KeyType>.push_back(std::make_pair(keys, keyNum));
            }

            // 为根节点分配键存储空间
            KeyType *rootKeys = NULL;
            my_alloc::align_zalloc((void **) &rootKeys, (seg_message.size() + 1) * sizeof(KeyType));
            keyNum = 0;

            // 收集并存储根节点的键信息
            rootKeys[keyNum++] = min_key;
            for (const auto &item: seg_message) {
                rootKeys[keyNum++] = item.key;
            }

            // 将根节点的键信息存储到treeLayout中
            treeLayout<KeyType>.push_back(std::make_pair(rootKeys, keyNum));
        }
/**
 根据索引布局自顶向下构建索引树。
 treelayout 索引树的布局信息，包含各层节点的关键字列表及节点数量。
 */
        void buildTree(TreeLayoutType &treelayout, const std::vector<KV> &kvs, innerNode<KeyType> *&root, int rootSize) {
            // 初始化最小键值、最大键值和树的高度。
            KeyType min_key = kvs[0].first;
            KeyType max_key = kvs.back().first;
            int treeHeight = treelayout.size();
            my_alloc::statistic.maxdegree=treeHeight;
            // 获取树的分支数。
            int fanout = treelayout.back().second;
            // 分配内存给根节点并初始化。
            root = reinterpret_cast<innerNode<KeyType> *>(Alloc::PMZAlloc(rootSize));
            root->setInnerNode(min_key, max_key, fanout);
            Alloc::Persist(&root, sizeof(innerNode<KeyType>));

            // 初始化父节点指针、范围起始点数组和子节点指针。
            double *parents_range_froms;
            my_alloc::allocMem(parents_range_froms, 2);
            parents_range_froms[0] =0;
            parents_range_froms[1] = max_key_;
            uintptr_t *parents;
            my_alloc::allocMem(parents, 1);
            uintptr_t *children = NULL;
            parents[0] = reinterpret_cast<uintptr_t> (root);
            int n_parents = 1;
            int act_total_N_children = 0;

            // 从倒数第二层开始向上构建索引的各个节点。
            for (int h = treeHeight - 1; h > 0; h--) {
                std::pair<KeyType *, int> &keys_list = treelayout[h - 1];
                auto _pair = create_children(h, parents, n_parents, parents_range_froms, keys_list,
                                             act_total_N_children);
                children = _pair.first;
                double *children_range_froms = _pair.second;
#ifndef LARGE_DRAM
                free(parents);
#endif
                parents = children;
#ifndef LARGE_DRAM
                free(parents_range_froms);
#endif
                parents_range_froms = children_range_froms;
                n_parents = act_total_N_children;
            }

            // 对所有叶子节点进行初始化。
            KeyType *all_keys = treelayout[0].first;
            int num = treelayout[0].second;
            for (long i = 0; i < num; ++i) {
                auto leaf = reinterpret_cast<std::pair<int, innerSlot<KeyType>*> *>(find_leaf(all_keys[i]));
                leaf->first++;
            }
            for (auto it: treelayout) {
#ifndef LARGE_DRAM
                free(it.first);
#endif
            }

            // 初始化叶节点结构体用于存储键值对。
            long start_idx = 0;
            int kvsPerBlock = (BLOCK_SIZE - sizeof(BlockHeader))/sizeof(KV);
            SegType *pre = NULL;
            my_alloc::statistic.segNUM = act_total_N_children;
            for (int i = 0; i < act_total_N_children; ++i) {
                auto pair=reinterpret_cast<std::pair<int, innerSlot<KeyType>*> *>(children[i]);
                uint32_t blockNum = pair->first/INIT_RATIO/kvsPerBlock+1;
                auto seg=reinterpret_cast<SegType *>(Alloc::PMZAlloc(sizeof(SegType)+blockNum*sizeof(Block<KeyType,ValueType>)));
                if(pre==NULL){
                    segments_head_=seg;
                    pre=seg;
                }
                else{
                    pre->set_next_segment(seg);
                    seg->set_pre_segment(pre);
                    pre=seg;
                }
                seg->init(blockNum, pair->first);
                seg->bulk_loading(kvs, start_idx,pair->second,version);
                Alloc::Persist(seg, sizeof(SegType)+blockNum*sizeof(Block<KeyType,ValueType>));
                start_idx += pair->first;
                delete pair;
            }
#ifndef LARGE_DRAM
            free(children);
    free(parents_range_froms);
#endif
        }
        //查找给定键值对应的叶子节点
        uint64_t find_leaf(KeyType key) { 
            innerNode<KeyType> *node = root;
            int idx = LR_PRED( root->slope,root->minkey, key, root->fanout);
            auto slot = &root->slots[idx];
            while (!slot->isleaf()) {
                node = reinterpret_cast<innerNode<KeyType> *>(slot->getChildNodePtr());
                idx = LR_PRED(slot->slope,slot->minkey, key, slot->size);
                slot = &node->slots[idx];
            }
            return slot->getChildNodePtr();
        }
        //查找给定键值对应的叶子节点
        SegType *find_seg(const KeyType &key, innerSlot<KeyType> * &_slot) {
            // 从根节点开始查找
            innerNode<KeyType> *node = root;
            int idx = LR_PRED(root->slope,root->minkey, key, root->fanout);
            auto slot = &root->slots[idx];
#ifndef LEVELNUM

            while (!slot->isleaf()) {
                node = reinterpret_cast<innerNode<KeyType> *>(slot->getChildNodePtr());
                idx = LR_PRED(slot->slope,slot->minkey, key, slot->size);
                slot = &node->slots[idx];
            }
#else
            int height=1;
           while (!slot->isleaf()) {
                height++;
                node = reinterpret_cast<innerNode<KeyType> *>(slot->getChildNodePtr());
                idx = LR_PRED(slot->slope,slot->minkey, key, slot->size);
                slot = &node->slots[idx];
            }
           std::cout<<height<<std::endl;
           //exit(-1);
#endif
            _slot=slot;
            return reinterpret_cast<SegType * >(slot->getChildNodePtr());
        }

/**
 * 创建子节点并初始化其范围和指向。
 * height 树的高度
 * parents 父节点的指针数组。
 * n_parents 父节点数组的大小。
 * parents_range_froms 父节点范围的起始值数组。
 * keys_list 包含键列表和键的数量。
 * act_total_N_children 实际子节点总数，输出参数。
 * 一个包含子节点指针和子节点范围起始值数组的std::pair。
 */

        std::pair<uintptr_t *, double *> create_children(const int &height, uintptr_t *parents, int n_parents,
                                                         double *parents_range_froms, std::pair<KeyType *, int> &keys_list,
                                                         int &act_total_N_children) {
            act_total_N_children = 0; // 初始化实际子节点总数为0

            // 计算实际子节点总数
            for (int i = 0; i < n_parents; ++i) {
                auto p = reinterpret_cast<innerNode<KeyType> *>(parents[i]);
                act_total_N_children += p->fanout;
            }

            uintptr_t *_children; // 子节点指针数组
            my_alloc::allocMem(_children, act_total_N_children); // 分配内存给子节点指针数组

            double *children_range_froms; // 子节点范围起始值数组
            if (height > 1) {
                my_alloc::allocMem(children_range_froms, act_total_N_children + 1); // 高度大于1时，分配内存给子节点范围数组
            }

            KeyType *list = keys_list.first; // 键列表
            int keys_num = keys_list.second; // 键列表的数量
            int last_idx = 0; // 上一个键的索引
            int cursor = 0; // 当前子节点的索引

            // 遍历父节点，为每个父节点创建子节点
            for (int i = 0; i < n_parents; ++i) {
                auto parent = reinterpret_cast<innerNode<KeyType> *>(parents[i]);
                int fanout = parent->fanout; // 父节点的子节点数量
                double range_from = parents_range_froms[i]; // 父节点范围的起始值
                double parent_range_to = parents_range_froms[i + 1]; // 父节点范围的结束值

                // 为每个子节点创建和初始化
                for (int child_id = 0; child_id < fanout; ++child_id) {
                    double range_to = 1.0 * (child_id + 1) / parent->slope + parent->minkey; // 计算子节点范围的结束值
                    // 调整范围结束值以确保不超过父节点范围
                    if (range_to > parent_range_to) {
                        assert(child_id >= (fanout - 1));
                        range_to = parent_range_to;
                    }
                    if (child_id == fanout - 1) {
                        range_to = parent_range_to;
                    }

                    int idx = array_lower_bound(list, range_to, 0, keys_num); // 查找子节点键的插入位置
                    int n_keys_this_child = idx - last_idx; // 计算当前子节点包含的键数量

                    if (height > 1) {
                        // 为内部节点分配内存并初始化
                        int nodefanout = n_keys_this_child >= 1 ? n_keys_this_child : 1;
                        auto child = reinterpret_cast<innerNode<KeyType> *>(Alloc::PMZAlloc(
                                sizeof(innerNode<KeyType>) + nodefanout * sizeof(innerSlot<KeyType>)));
                        child->setInnerNode(list + last_idx, n_keys_this_child);
                        Alloc::Persist(child, sizeof(innerNode<KeyType>));
                        children_range_froms[cursor] = range_from;
                        auto ptr = reinterpret_cast<uintptr_t>(child);
                        _children[cursor++] = ptr;
                        parent->slots[child_id].setSlot(child->minkey, ptr, child->slope, child->fanout, 0, 0);
                        Alloc::Persist(&parent->slots[child_id], sizeof(innerSlot<KeyType>));
                    } else {
                        // 为叶节点分配内存并初始化
                        auto seg = new std::pair<int, innerSlot<KeyType>*>(0, &parent->slots[child_id]);
                        auto ptr = reinterpret_cast<uintptr_t>(seg);
                        _children[cursor++] = ptr;
                        parent->slots[child_id].setSlot(ptr, version, leafset);
                    }

                    last_idx = idx; // 更新上一个键的索引
                    range_from = range_to; // 更新范围的起始值
                }
            }

            if (height > 1) {
                children_range_froms[cursor] = parents_range_froms[n_parents]; // 更新子节点范围数组的最后一个值
            }

            return std::make_pair(_children, children_range_froms); // 返回子节点指针和范围数组
        }


        //自底向上构建索引布局和自顶向下构建索引        // Keys must be sorted.
        void BulkLoad(const std::vector<KV> &kvs) {
            min_key_ = std::min(min_key_, kvs.front().first);
            max_key_ = std::max(max_key_, kvs.back().first);
            setTreeLayout(kvs);
            buildTree(treeLayout<KeyType>, kvs, root, 256);
        }
        //根据误差error自底向上构建索引布局和自顶向下构建索引        // Keys must be sorted.
        void BulkLoad(const std::vector<KV> &kvs,size_t error) {
            min_key_ = std::min(min_key_, kvs.front().first);
            max_key_ = std::max(max_key_, kvs.back().first);
            setTreeLayout(kvs,error);
            buildTree(treeLayout<KeyType>, kvs, root, 256);
        }
        //索引恢复检查操作
        void recoveryCheck(){
            auto seg=segments_head_;
            Segment<KeyType, ValueType> * pre=NULL;
            while(seg!=NULL){
                seg->set_pre_segment(pre);
                pre=seg;
                auto _slot=seg->getSlot();
                blockMeta *metas=NULL;
                my_alloc::allocMem(metas, seg->getNum());
                seg->set_metas(metas);
                _slot->v=version;
                _slot->setBufferPtr(reinterpret_cast<uint64_t>(metas));
                _slot->release_lock();
                if(seg->getSlot()->getChildNodePtr()!=reinterpret_cast<uint64_t>((seg))){
                    _slot->setChildNodePtr(reinterpret_cast<uintptr_t>(seg));
                    //slot->setBufferPtr(reinterpret_cast<uint64_t>(metas));
                    _slot->size = seg->getNum();
                    _slot->minkey = seg->get_minKey();
                    _slot->slope = seg->get_slope();
                    Alloc::Persist(_slot, sizeof(innerSlot<KeyType>));
                }
                seg=seg->get_next_segment();
            }
        }
        void freeSegMetas(){
            auto seg=segments_head_;
            while(seg!=NULL){
                seg->freeMetas();
                seg=seg->get_next_segment();
            }

        }


/**
 插入一个新的键值对，如果键已存在，则更新其值。
 如果插入或更新成功，则返回 true；如果插入失败，则返回 false。
 */
        bool insert(const KeyType &key, ValueType value) {
            innerSlot<KeyType> *slot=NULL; // 初始化指向叶节点的槽位的指针为 NULL
            RETRY:
            auto seg=find_seg(key,slot); // 查找键对应的叶节点
            auto ret = seg->upsert(key, value,slot); // 尝试插入或更新操作，并获取操作结果

            // 判断操作结果并处理
            if (ret == OpStatus::Inserted||ret==OpStatus::Updated) {
                // 如果键已插入或已更新，则成功
                return true;
            };

            if (ret == OpStatus::Failed) {
                // 如果操作失败，则直接返回 false
                return false;
            };
            if(ret==OpStatus::lock){
                // 如果操作因锁定而失败，则重试
                goto RETRY;
            }

            if (ret == OpStatus::SMO) {
                // 如果操作需要结构化移动优化（SMO），则重对叶节点进行扩展
                Retrain(seg,slot);
                goto RETRY;
            }
        }
        /**
        在索引中搜索指定的键，并返回value
        如果找到了键，返回true；如果未找到，返回false。
       */
        bool search(const KeyType &key, ValueType &value) {
            innerSlot<KeyType> *slot=NULL;
            RETRY:
#ifndef INNERTIME
            //根据可以查找键
            auto ret= find_seg(key,slot)->find(key, value,slot);
            if(ret==OpStatus::lock||ret==OpStatus::Retry){
                goto RETRY;// 如果操作被锁定或需要重试，则重试查找
            }
            if(ret==OpStatus::Find){
                return true; // 如果成功找到键，返回true。
            }
            if(ret==OpStatus::Failed||ret==OpStatus::NotFound){
                return false; // 如果查找失败或键不存在，返回false
            }
#else//只测试内部节点查找延迟
            find_seg(key,slot);
            return true;
#endif
        }
//键值对的删除用插入DELETE标记的方式进行
        bool erase(const KeyType &key) {
            innerSlot<KeyType> *slot=NULL;
            RETRY:
            auto seg=find_seg(key,slot);
            auto ret =seg->erase(key,slot);
            if (ret == OpStatus::Inserted||ret==OpStatus::Updated) {
                return true;
            };
            if (ret == OpStatus::Failed) {
                return false;
            };
            if(ret==OpStatus::lock){
                goto RETRY;
            }
            if (ret == OpStatus::SMO) {
                Retrain(seg,slot);
                goto RETRY;
            }
        }

        bool update(KeyType key, ValueType value) {
             return false;
        }
/**
执行范围查找操作。
该函数用于在数据结构中进行范围查找，从给定的起始键开始，扫描指定数量的元素，并将结果存储在一个map中。
这个操作会跨越多个叶节点，如果起始键位于一个叶节点的中间，它会从那个数据块开始扫描。
start_key 扫描的起始键。
size 扫描希望获取的元素数量。
results 用于存储扫描结果的map，键为KeyType，值为ValueType。
return 实际扫描到的元素数量。
 */
        int rangeScan(const KeyType &start_key, const int &size, std::map<KeyType, ValueType> &results) {
            int idx=-1;
            int blockNum=0;
            int scannum = 0;
            int deleted = 0;//已标记为删除的元素数量
            // 找到起始键所在的叶节点，并开始范围查找
            auto seg=find_seg(start_key,idx,blockNum);
            scannum = seg->rangeScan(start_key, size, deleted, results, idx,blockNum);
            // 如果已扫描到的元素数量小于期望的大小，则继续扫描下一个叶节点
            while (scannum < size && (seg = seg->get_next_segment())) {
                // 在每个叶节点的开始处进行扫描，累加扫描到的元素数量
                scannum += seg->scanAtBegin(start_key, size, deleted, results, seg->getBlockNum());
            }
            return scannum;
        }


        SegType *GetSegmentHead() {
            return segments_head_;
        }

        SegType *GetSegmentTail() {
            return segments_tail_;
        }

        void closePool() {
            // Alloc::_pop.ClosePool();
        }
    private:
/**
  对给定的叶节点进行扩展操作。
  segment 指向当前需要重新训练的叶节点的指针。
  slot 指向当前叶节点中的内部插槽（innerSlot）的指针，用于存储指向的叶节点信息。
  通过创建一个新的叶节点来存储键值对
 */
        void Retrain(SegType *segment,innerSlot<KeyType>* _slot) {
            // 尝试获取叶节点扩展分割锁，如果获取失败，则直接返回
            if (!_slot->try_get_split_lock()) {
                return;
            }
            // 尝试获取当前叶节点的左兄弟叶节点锁，如果不存在或获取锁失败，则释放锁并返回
            SegType * left_sibling = segment->pre_segment();
            if((!left_sibling)||(!left_sibling->getSlot()->try_get_split_lock())){
                _slot->release_lock();
                return;
            }

            // 尝试获取当前叶节点的右兄弟叶节点锁，如果不存在或获取锁失败，则释放已获取的锁并返回
            auto right_sibling = segment->get_next_segment();
            if((!right_sibling)||(!right_sibling->getSlot()->try_get_split_lock())){
                _slot->release_lock();
                if(left_sibling)
                    left_sibling->getSlot()->release_lock();
                return;
            }

            // 成功获取了所有相关叶节点锁的锁
            _slot->get_write_lock();

            // 计算新的叶节点锁中块（Block）的数量
            int blockSize = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
            std::map<KeyType, ValueType> results;
            int kv_num = segment->getAllKVs(results,_slot->size);
            uint32_t blockNum = kv_num/INIT_RATIO/blockSize+1;

            // 分配内存给新的叶节点，并初始化
            auto new_seg=reinterpret_cast<SegType *>(Alloc::PMZAlloc(sizeof(SegType)+blockNum*sizeof(Block<KeyType,ValueType>)));
            new_seg->init(blockNum,kv_num);
            new_seg->resize(results,_slot,segment->pre_segment(),segment->get_next_segment());

            // 持久化新的叶节点
            Alloc::Persist(new_seg, sizeof(SegType)+blockNum*sizeof(Block<KeyType,ValueType>));

            // 设置读锁，并更新插槽与叶节点的相关信息
            _slot->get_read_lock();
            if(right_sibling)
                right_sibling->set_pre_segment(new_seg);
            if(left_sibling)
                left_sibling->switch_next_segment(new_seg);

            // 更新插槽大小、最小键值和斜率
            _slot->size = blockNum;
            _slot->minkey = new_seg->get_minKey();
            _slot->slope = new_seg->get_slope();
            _slot->setChildNodePtr(reinterpret_cast<uintptr_t>(new_seg));

            // 持久化插槽信息
            Alloc::Persist(_slot, sizeof(innerSlot<KeyType}));

        // 释放原叶节点的元数据，并释放相关叶节点的锁
        segment->freeMetas();
        if(right_sibling)
        right_sibling->getSlot()->release_lock();
        if(left_sibling)
        left_sibling->getSlot()->release_lock();
        _slot->release_lock();
    }
        KeyType min_key_;
        KeyType max_key_;
        size_t max_error_;
        SegType *segments_head_, *segments_tail_;
        innerNode<KeyType> *root;
        uint8_t version;
        bool crash;
    };

} 
#endif //PLTREE_H
