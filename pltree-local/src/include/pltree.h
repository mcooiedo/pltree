//
// Created by mcooiedo on 22-12-4.
//

#ifndef PLTREE_H
#define PLTREE__H

#include <algorithm>
#include <cmath>
#include <iostream>
#include "segment.h"
#include "allocator.h"
//#include <mcheck.h>
//#define BTREE
#define ROOT_SIZE 7
#define  LEVEL_SIZE 20

namespace plt {

 template<typename KeyType>
     class RMI{
  
    public :
        // static const int ROOT_SIZE=8;
        int rootsize;
        size_t max_error;
        innerSlot<KeyType> root[ROOT_SIZE];
        RMI(size_t max_error =8):max_error(max_error),rootsize(0){

        }
        void buildTree(std::vector<innerSlot<KeyType>> &sis){
            std::vector<innerSlot<KeyType>> tmp=sis;
            std::vector<innerSlot<KeyType>> newtmp;
//uint8_t level = 1;
            while ( tmp.size() > ROOT_SIZE) {
                plt::Builder<KeyType> asb(sis.front().maxkey, sis.back().maxkey, max_error);
                for (const auto& item : tmp) {
                    asb.AddKey(item.maxkey);
                }
                asb.Finalize();
                auto &seg_message = asb.get_segments_message();
                for (const SegmentMessage<KeyType> & msg : seg_message) {
                     auto *node=reinterpret_cast<innerNode<KeyType> *>(Alloc::PMZAlloc(msg.size*sizeof(innerSlot<KeyType>)+sizeof(innerNode<KeyType>)));
                    node->addKV(msg,tmp);
                     Alloc::Persist(node, msg.size*sizeof(innerSlot<KeyType>)+sizeof(innerNode<KeyType>));
                    //node->level=level;
                    innerSlot<KeyType>  slot;
                    slot.setSlot(tmp[msg.offset].maxkey,msg.key,(uintptr_t)node,msg.slope,msg.size,0,0);
                    newtmp.push_back(slot);
                }
               // level++;
                tmp=newtmp;
                newtmp.clear();
            }
            for(int i=0;i<tmp.size();++i){
                root[i]=tmp[i];
                rootsize++;
            }
            //root=tmp;
        }
        KeyType getMaxKey(){
            root.back().first;
        }
        innerSlot<KeyType>* find(KeyType key){
            int i=0;
            for(;i<rootsize;++i){
                if(root[i].maxkey>=key)
                    break;
            }
            auto slot=&root[i];
           // auto segInf=&root[i];
            while (!slot->isleaf()){
                auto node = reinterpret_cast<innerNode<KeyType> *>(slot->getChildNodePtr());
                slot=reinterpret_cast<innerSlot<KeyType>*>(node->find(key,slot,max_error));
            }
            return slot;
        }
         uint64_t  find_seg(KeyType key, innerSlot<KeyType>* &_slot){
             int i=0;
             for(;i<rootsize;++i){
                 if(root[i].maxkey>=key)
                     break;
             }
             auto slot=&root[i];
             // auto segInf=&root[i];
             while (!slot->isleaf()){
                 auto node = reinterpret_cast<innerNode<KeyType> *>(slot->getChildNodePtr());
                 slot=reinterpret_cast<innerSlot<KeyType>*>(node->find(key,slot,max_error));
             }
                _slot=slot;
             return slot->getChildNodePtr();
         }

    };

    template<typename KeyType, typename ValueType, typename Alloc=my_alloc::MemPoolAllocator>//PMAllocator>
    class PLTree {
        typedef struct KV<KeyType, ValueType> KV;
        //typedef  BufferNode<KeyType,ValueType> BufferNodeType;
        typedef Segment<KeyType, ValueType> SegType;
        typedef std::vector<std::pair<KeyType *, int>> TreeLayoutType;

    public:
        PLTree(bool checkrecover=true)
               {
            //tree = RMI<KeyType>(8);
            if (Alloc::isRecover()&&checkrecover) {
                version=version+1==0?1:version+1;
                recoveryCheck();
            } else {
                version=1;
                min_key_=std::numeric_limits<KeyType>::max();
                max_key_=std::numeric_limits<KeyType>::min();
                max_error_=256;
                segments_head_=nullptr;
                segments_tail_=nullptr;
                //global_buffer=new(Alloc::PMZAlloc(sizeof(BufferNodeType))) BufferNodeType();
                //Alloc::Persist(&global_buffer,sizeof (BufferNodeType*));
                //# SegType::segment_allocated_byte = 0;
                //global_overflow_buffer_= new OverflowBuffer<KeyType, ValueType>;//###
                //### global_overflow_buffer_= new(Alloc::pmAlloc(sizeof(OverflowBuffer<KeyType, ValueType>))) OverflowBuffer<KeyType, ValueType>(alloc);
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


        void setTreeLayout(const std::vector<KV> &kvs) {
            KeyType min_key = kvs[0].first;
            KeyType max_key = kvs.back().first;
            plt::Builder<KeyType> asb(min_key, max_key, max_error_/*,alloc*/);
            KeyType *keys = NULL;
            my_alloc::align_zalloc((void **) &keys, kvs.size() * sizeof(KeyType));
            int keyNum = 0;
            for (const auto &kv: kvs) {
                asb.AddKey(kv.first);
                keys[keyNum++] = kv.first;
            }
            asb.Finalize();
            auto &seg_message = asb.get_segments_message();
            treeLayout<KeyType>.push_back(std::make_pair(keys, keyNum));
            int error=256;
            while (seg_message.size() + 1 > ROOT_SIZE) {
                //error=error>256?256:error*8;
                KeyType *keys = NULL;
                my_alloc::align_zalloc((void **) &keys, (seg_message.size() + 1) * sizeof(KeyType));
                int keyNum = 0;
                plt::Builder<KeyType> asb(min_key, seg_message.back().key, error);
                asb.AddKey(min_key);
                keys[keyNum++] = min_key;
                for (const auto &item: seg_message) {
                    asb.AddKey(item.key);
                    keys[keyNum++] = item.key;
                }
                asb.Finalize();

                seg_message = asb.get_segments_message();
                treeLayout<KeyType>.push_back(std::make_pair(keys, keyNum));
            }
            KeyType *rootKeys = NULL;//=new KeyType[seg_message.size()+1];
            my_alloc::align_zalloc((void **) &rootKeys, (seg_message.size() + 1) * sizeof(KeyType));
            keyNum = 0;
            rootKeys[keyNum++] = min_key;
            for (const auto &item: seg_message) {
                rootKeys[keyNum++] = item.key;
            }
            treeLayout<KeyType>.push_back(std::make_pair(rootKeys, keyNum));
        }

        void buildTree(TreeLayoutType &treelayout, const std::vector<KV> &kvs, innerNode<KeyType> *&root, int rootSize) {
            KeyType min_key = kvs[0].first;
            KeyType max_key = kvs.back().first;
            int treeHeight = treelayout.size();
            my_alloc::statistic.maxdegree=treeHeight;
            int fanout = treelayout.back().second;
            root = reinterpret_cast<innerNode<KeyType> *>(Alloc::PMZAlloc(rootSize));
            root->setInnerNode(min_key, max_key, fanout);
            Alloc::Persist(&root, sizeof(innerNode<KeyType>));
            double *parents_range_froms;
            my_alloc::allocMem(parents_range_froms, 2);
            parents_range_froms[0] =0;//  min_key_;
            parents_range_froms[1] = max_key_;
            //innerNode<KeyType> **parents = new innerNode<KeyType>*[1];
            //innerNode<KeyType> **children = NULL;
            uintptr_t *parents;// = new uintptr_t[1];
            my_alloc::allocMem(parents, 1);
            //my_alloc::align_zalloc((void**) &parents,sizeof(uintptr_t));
            uintptr_t *children = NULL;
            parents[0] = reinterpret_cast<uintptr_t> (root);
            int n_parents = 1;
            int act_total_N_children = 0;
            for (int h = treeHeight - 1; h > 0; h--) {
                std::pair<KeyType *, int> &keys_list = treelayout[h - 1];
                auto _pair = create_children(h, parents, n_parents, parents_range_froms, keys_list,
                                             act_total_N_children);
                children = _pair.first;
                //auto p= reinterpret_cast<innerNode<KeyType>*>(children[0]);
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
            KeyType *all_keys = treelayout[0].first;
            int num = treelayout[0].second;
            for (long i = 0; i < num; ++i) {
/*                if(all_keys[i]==5067167772873796717){
                    int y=0;
                }*/
                auto leaf = reinterpret_cast<std::pair<int, innerSlot<KeyType>*> *>(find_leaf(all_keys[i]));
                leaf->first++;
            }
            for (auto it: treelayout) {
#ifndef LARGE_DRAM
                free(it.first);
#endif
            }
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
                //cout<<i<<endl;
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

        uint64_t find_leaf(KeyType key) { ///??
            //auto slot = root->findSlot(key);
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

        innerSlot<KeyType> *find_slot(KeyType key, int &slotChildNum) {
            innerNode<KeyType> *node = root;
            int idx = LR_PRED(root->slope,root->minkey, key, root->fanout);
            auto slot = &root->slots[idx];
            while (!slot->isleaf()) {
                node = reinterpret_cast<innerNode<KeyType> *>(slot->getChildNodePtr());
                idx = LR_PRED(slot->slope,slot->minkey, key, slot->size);
                slot = &node->slots[idx];
            }
            /*while (slot->childNum < 0 && idx >= 0) {
                idx--;
                slot = &node->slots[idx];
            }*/
            slotChildNum = slot->childNum;
            return slot;
        }
        SegType *find_seg(const KeyType &key, int &blockIdx,blockMeta * &ptr,innerSlot<KeyType> * &_slot) {
            innerNode<KeyType> *node = root;
            int idx = LR_PRED(root->slope,root->minkey, key, root->fanout);
            auto slot = &root->slots[idx];
            while (!slot->isleaf()) {
                node = reinterpret_cast<innerNode<KeyType> *>(slot->getChildNodePtr());
                idx = LR_PRED(slot->slope,slot->minkey, key, slot->size);
                slot = &node->slots[idx];
            }
            ptr=reinterpret_cast<blockMeta * >(slot->getBufferPtr());
            _slot=slot;
            auto block=(key-slot->minkey) * slot->slope /15;
            //slot->predict_block(key,blockIdx,15);
            blockIdx=block<0?0:(block>=slot->size?slot->size-1:block);
            return reinterpret_cast<SegType * >(slot->getChildNodePtr());
        }
        SegType *find_seg(const KeyType &key, int &blockIdx,blockMeta * &ptr) {
            innerNode<KeyType> *node = root;
            int idx = LR_PRED(root->slope,root->minkey, key, root->fanout);
            auto slot = &root->slots[idx];
            while (!slot->isleaf()) {
                node = reinterpret_cast<innerNode<KeyType> *>(slot->getChildNodePtr());
                idx = LR_PRED(slot->slope,slot->minkey, key, slot->size);
                slot = &node->slots[idx];
            }
            ptr=reinterpret_cast<blockMeta * >(slot->getBufferPtr());
            auto block=(key-slot->minkey) * slot->slope /15;
            //slot->predict_block(key,blockIdx,15);
            blockIdx=block<0?0:(block>=slot->size?slot->size-1:block);
            return reinterpret_cast<SegType * >(slot->getChildNodePtr());
        }
        SegType *find_seg(const KeyType &key, int &blockIdx,int &blockNum) {
            innerNode<KeyType> *node = root;
            int idx = LR_PRED(root->slope,root->minkey, key, root->fanout);
            auto slot = &root->slots[idx];
            while (!slot->isleaf()) {
                node = reinterpret_cast<innerNode<KeyType> *>(slot->getChildNodePtr());
                idx = LR_PRED(slot->slope,slot->minkey, key, slot->size);
                slot = &node->slots[idx];
            }
            blockNum=slot->size;
            auto block=(key-slot->minkey) * slot->slope /15;
            //slot->predict_block(key,blockIdx,15);
            blockIdx=block<0?0:(block>=slot->size?slot->size-1:block);
            return reinterpret_cast<SegType * >(slot->getChildNodePtr());
        }


        std::pair<uintptr_t *, double *> create_children(const int &height, uintptr_t *parents, int n_parents,
                                                    double *parents_range_froms, std::pair<KeyType *, int> &keys_list,
                                                    int &act_total_N_children) {
            act_total_N_children = 0;
            for (int i = 0; i < n_parents; ++i) {
                /*    if (parents[i] == 0)
                         continue;*/
                auto p = reinterpret_cast<innerNode<KeyType> *>(parents[i]);
                act_total_N_children += p->fanout;
            }
            uintptr_t *_children;// = new uintptr_t[act_total_N_children];
            my_alloc::allocMem(_children, act_total_N_children);

            double *children_range_froms;
            if (height > 1) {
                my_alloc::allocMem(children_range_froms, act_total_N_children + 1);
                //children_range_froms = new double[act_total_N_children];
            }
            KeyType *list = keys_list.first;
            int keys_num = keys_list.second;

            int last_idx = 0;
            int cursor = 0;
            for (int i = 0; i < n_parents; ++i) {
                /*        if (parents[i] == 0)
                             continue;*/
                auto parent = reinterpret_cast<innerNode<KeyType> *>(parents[i]);
                //innerNode<KeyType> *parent_1 = reinterpret_cast<innerNode<KeyType> *>(parents[i - 1]);
                int fanout = parent->fanout;
                double range_from = parents_range_froms[i];
                double parent_range_to = parents_range_froms[i + 1];
                for (int child_id = 0; child_id < fanout; ++child_id) {
/*                    if(parent->slope==0)
                        int y=0;*/
                    double range_to = 1.0 * (child_id + 1) / parent->slope+parent->minkey ;
                    if (range_to > parent_range_to) {
                        assert(child_id >= (fanout - 1));
                        range_to = parent_range_to;
                    }
                    if (child_id == fanout - 1) {
                        range_to = parent_range_to;
                    }

                    int idx = array_lower_bound(list, range_to, 0, keys_num);
                    int n_keys_this_child = idx - last_idx;
                    if (height > 1) {
                        int nodefanout=n_keys_this_child>=1?n_keys_this_child:1;
                        auto child = reinterpret_cast<innerNode<KeyType> *>(Alloc::PMZAlloc(
                                sizeof(innerNode<KeyType>) +nodefanout * sizeof(innerSlot<KeyType>)));
                        child->setInnerNode(list + last_idx, n_keys_this_child);
                        Alloc::Persist(child, sizeof(innerNode<KeyType>));
                        children_range_froms[cursor] = range_from;
                        auto ptr = reinterpret_cast<uintptr_t>(child);
                        //auto p=reinterpret_cast<innerNode<KeyType> *>(ptr);
                        _children[cursor++] = ptr;
                        parent->slots[child_id].setSlot(child->minkey, ptr, child->slope, child->fanout, 0, 0);
                        Alloc::Persist(&parent->slots[child_id], sizeof(innerSlot<KeyType>));
                        /*                      if (n_keys_this_child >= 1) {
                                                  auto child = reinterpret_cast<innerNode<KeyType> *>(Alloc::PMZAlloc(
                                                          sizeof(innerNode<KeyType>) + n_keys_this_child * sizeof(innerSlot<KeyType>)));
                                                  child->setInnerNode(list + last_idx, n_keys_this_child);
                                                  children_range_froms[cursor] = range_from;
                                                  auto ptr = reinterpret_cast<uintptr_t>(child);
                                                  //auto p=reinterpret_cast<innerNode<KeyType> *>(ptr);
                                                  _children[cursor++] = ptr;
                                                  parent->slots[child_id].setSlot(child->b, ptr, child->slope, child->fanout, 1, 0);
                                              } else {

                            auto child = reinterpret_cast<innerNode<KeyType> *>(Alloc::PMZAlloc(
                                    sizeof(innerNode<KeyType>) + sizeof(innerSlot<KeyType>)));
                            child->setInnerNode(list + last_idx, n_keys_this_child);
                            children_range_froms[cursor] = range_from;
                            auto ptr = reinterpret_cast<uintptr_t>(child);
                            _children[cursor++] = ptr;
                            parent->slots[child_id].setSlot(child->b, ptr, child->slope, child->fanout, 1, 0); */
                            /*  children_range_froms[cursor] = range_from;
                              _children[cursor++] = 0;
                              if (child_id == 0) {
                                  int t = cursor - 1;
                                  while (_children[t] == 0) {
                                      t--;
                                  }
                                  auto p = reinterpret_cast<innerNode<KeyType> *>(_children[t]);
                                  parent->slots[child_id].setSlot(p->minkey, _children[t],
                                                                  p->slope, p->fanout, -1,
                                                                  0);
                              } else
                                  parent->slots[child_id].setSlot(parent->slots[child_id - 1].minkey,
                                                                  parent->slots[child_id - 1].childNodePtr,
                                                                  parent->slots[child_id - 1].slope,
                                                                  parent->slots[child_id - 1].size, -1,
                                                                  0);
                        }*/
                    } else {
                        auto seg =  new std::pair<int, innerSlot<KeyType>*>(0,&parent->slots[child_id]);//reinterpret_cast<SegType *>(Alloc::PMZAlloc(sizeof(SegType)));
                        //seg->setSlot(&parent->slots[child_id]);
                        auto ptr = reinterpret_cast<uintptr_t>(seg);
                        _children[cursor++] = ptr;
                        parent->slots[child_id].setSlot(ptr, version, leafset);
                        /*if (n_keys_this_child >= 1) {
                             auto seg = reinterpret_cast<SegType *>(Alloc::PMZAlloc(sizeof(SegType)));
                             seg->setSlot(&parent->slots[child_id]);
                             auto ptr = reinterpret_cast<uintptr_t>(seg);
                             _children[cursor++] = ptr;
                             if(ptr==140692395202496)
                                 int y=0;
                             parent->slots[child_id].setSlot(ptr, 1, 1);
                         } else {
                             _children[cursor++] = 0;
                             if (child_id == 0) {
                                 int t = cursor - 2;
                                 while (_children[t] == 0) {
                                     t--;
                                 }
                                 parent->slots[child_id].setSlot(_children[t], -1,1);
                                 if(_children[t]==140692395202496)
                                     int y=0;
                             } else{
                                 parent->slots[child_id].setSlot(_children[cursor-2], -1, 1);
                                 if(_children[cursor-2]==140692395202496)
                                     int y=0;
                             }
                         }*/
                    }
                    last_idx = idx;
                    range_from = range_to;
                }
            }
            if (height > 1) {
                children_range_froms[cursor] = parents_range_froms[n_parents];
            }

            return std::make_pair(_children, children_range_froms);
        }


        // Keys must be sorted.
       /* void BulkLoad(const std::vector<KV> &kvs) {
            min_key_ = std::min(min_key_, kvs.front().first);
            max_key_ = std::max(max_key_, kvs.back().first);
            setTreeLayout(kvs);
            buildTree(treeLayout<KeyType>, kvs, root, 256);
        }*/
        void BulkLoad(const std::vector<KV> &kvs) {
            assert(kvs.size() > 0);
            //assert(keys.size() == values.size());
            // Build Alert.
            min_key_ = std::min(min_key_, kvs.front().first);
            max_key_ = std::max(max_key_, kvs.back().first);

            plt::Builder<KeyType> asb(min_key_, max_key_, max_error_/*,alloc*/);
            for (const auto& kv : kvs) {
                asb.AddKey(kv.first);
            }
            asb.Finalize();

            auto &seg_message = asb.get_segments_message();
            my_alloc::statistic.segNUM=seg_message.size();
            SegType *pre_seg = nullptr;
            std::vector<innerSlot<KeyType>> sis;//tree
            for (const SegmentMessage<KeyType> & msg : seg_message) {
                auto ptr=Alloc::PMZAlloc(sizeof(SegType)+(msg.size/INIT_RATIO/15+2)*sizeof(Block<KeyType,ValueType>));
                auto seg=reinterpret_cast<SegType *>(ptr);
                innerSlot<KeyType> * seginf;
                seg->AddKV(msg, kvs,seginf,0);
                seg->set_pre_segment(pre_seg);
                //seg->segbulkPersisit();
                if (pre_seg) pre_seg->set_next_segment(seg);
                else {
                    segments_head_ = seg;
                    Alloc::Persist(segments_head_, sizeof(SegType *));
                }
                pre_seg = seg;
               /// tree_.Insert(msg.key, reinterpret_cast<uintptr_t>(segMeta));
                sis.push_back(*seginf);
                delete seginf;//SegmentInf<KeyType>{msg.key, reinterpret_cast<uintptr_t>(segMeta)}
                // tree_.insert(msg.key, reinterpret_cast<uintptr_t>(segMeta));
            }
            tree.buildTree(sis);
            segments_tail_ = pre_seg;
            Alloc::Persist(segments_tail_, sizeof(SegType *));
//            tree_.print_node_msg();
        }

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

         bool insert(const KeyType &key, ValueType value) {

            int idx=-1;
            blockMeta *bm;
            innerSlot<KeyType> *slot=NULL;
            RETRY:
            auto seg=find_seg(key,idx,bm,slot);
            auto ret = seg->upsert(key, value, idx,bm,slot);
            if (ret == OpStatus::Inserted||ret==OpStatus::Updated) {
                //num_total_keys_ += 1;
                return true;
            };

            if (ret == OpStatus::Failed) {
                //num_total_keys_ += 1;
                return false;
            };
            if(ret==OpStatus::lock){
                goto RETRY;
            }

            if (ret == OpStatus::SMO) {
                Retrain(seg,slot);
                goto RETRY;
               ///my_alloc::statistic.retry++;
            }
           // return true;
            // return 1;
        }

        bool search(const KeyType &key, ValueType &value) {

          //  int idx=-1;
          //  auto seg = reinterpret_cast<SegType *>();
            innerSlot<KeyType> *slot=NULL;
            //blockMeta *bm;
        
            RETRY:
          //  slot=tree.find(key);
            auto seg=reinterpret_cast<SegType *>(tree.find_seg(key,slot));
            if(seg==NULL)
                    return false;
            return true;

            auto  ret= seg->find(key, value, slot);
            
            if(ret==OpStatus::lock||ret==OpStatus::Retry){
                goto RETRY;
            }
            if(ret==OpStatus::Failed||ret==OpStatus::NotFound){
                return false;
            }
            if(ret==OpStatus::Find){
                return true;
            }

        }

        bool erase(const KeyType &key) {
            int idx=-1;
            blockMeta *bm;
            innerSlot<KeyType> *slot=NULL;
            RETRY:
            auto seg=find_seg(key,idx,bm,slot);
            auto ret =seg->erase(key, idx,bm,slot);
            if (ret == OpStatus::Inserted||ret==OpStatus::Updated) {
                //num_total_keys_ += 1;
                return true;
            };

            if (ret == OpStatus::Failed) {
                //num_total_keys_ += 1;
                return false;
            };
            if(ret==OpStatus::lock){
                goto RETRY;
            }

            if (ret == OpStatus::SMO) {
                Retrain(seg,slot);
                goto RETRY;
                ///my_alloc::statistic.retry++;
            }


        }

        bool update(KeyType key, ValueType value) {
            /*if (__glibc_unlikely(segments_head_ == nullptr || key > max_key_)) {
                // return global_buffer->updateInNode(key,value);
            }
            auto seginf = GetSplineSegment(key);
            auto segmeta = reinterpret_cast<SegType *>(seginf->node);
            auto ret = segmeta->upsert(key, value, seginf);
            if (ret == OpStatus::Updated) return true;
            else*/ return false;
        }

        int rangeScan(const KeyType &start_key, const int &size, std::map<KeyType, ValueType> &results) {
            int idx=-1;
            int blockNum=0;
            int scannum = 0;
            int deleted = 0;

            auto seg=find_seg(start_key,idx,blockNum);
            scannum = seg->rangeScan(start_key, size, deleted, results, idx,blockNum);

            while (scannum < size && (seg = seg->get_next_segment())) {
                scannum += seg->scanAtBegin(start_key, size, deleted, results, seg->getBlockNum());
            }
            return scannum;
        }

        // Returns the size in bytes.
        size_t GetSizeInByte() const {

            //  return sizeof(*this) +  tree_.size() + SegType::segment_allocated_byte;
            // return sizeof(*this) +  tree_.size() + SegType::segment_allocated_byte;
            return 0;

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

        void Retrain(SegType *segment,innerSlot<KeyType>* _slot) {
            if (!_slot->try_get_split_lock()) {
                //std::cout<<"release_lock:"<<last_accelerator<<std::endl;
                return;
            }
            SegType * left_sibling = segment->pre_segment();
            if((!left_sibling)||(!left_sibling->getSlot()->try_get_split_lock())){
                _slot->release_lock();
                return;
            }
            auto right_sibling = segment->get_next_segment();
            if((!right_sibling)||(!right_sibling->getSlot()->try_get_split_lock())){
                _slot->release_lock();
                if(left_sibling)
                    left_sibling->getSlot()->release_lock();
                return;
            }
            _slot->get_write_lock();
            int blockSize = (BLOCK_SIZE - sizeof(BlockHeader)) / sizeof(KV);
            std::map<KeyType, ValueType> results;
            int kv_num = segment->getAllKVs(results,_slot->size);
            uint32_t blockNum = kv_num/INIT_RATIO/blockSize+1;
            auto new_seg=reinterpret_cast<SegType *>(Alloc::PMZAlloc(sizeof(SegType)+blockNum*sizeof(Block<KeyType,ValueType>)));
            new_seg->init(blockNum,kv_num);
            new_seg->resize(results,_slot,segment->pre_segment(),segment->get_next_segment());
            Alloc::Persist(new_seg, sizeof(SegType)+blockNum*sizeof(Block<KeyType,ValueType>));
            _slot->get_read_lock();
            if(right_sibling)
                right_sibling->set_pre_segment(new_seg);
           if(left_sibling)
                left_sibling->switch_next_segment(new_seg);

            //log free seg  连接new_seg

            _slot->setChildNodePtr(reinterpret_cast<uintptr_t>(this));
            //slot->setBufferPtr(reinterpret_cast<uint64_t>(metas));
            _slot->size = blockNum;
            _slot->minkey = new_seg->get_minKey();
            _slot->slope = new_seg->get_slope();
            Alloc::Persist(_slot, sizeof(innerSlot<KeyType>));
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
        RMI<KeyType> tree;
    };

} 
#endif //PLTREE_H
