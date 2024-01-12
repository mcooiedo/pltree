

#ifndef PLTREE_TEST_COMMON_H
#define PLTREE_TEST_COMMON_H

#include <cstddef>
#include <cstdint>
#include "utils.h"
//#include "bucket.h"
#define INIT_ENTRY_BIT 1
#define FANOUT 1
#define BUCKET_PER_ENTRY 3
#define FREE 0
#define DELETED 0
#define ENABLE_LBM 1
#define THRESHOLD 50
#define THRESHOLD1 -1
#define BITS_PER_ITEM 12
#define LEVEL_BITS 4
#define FILTER_ITEMS 1000
#define FP_LEN 31
#define BUCKET_LEN 16

#define INIT_RATIO 0.5
#define SCALE_RATIO 0.5
#define BLOCK_SIZE 256


#define SIZE 15
#define HASH_SEED  0xDEADBEEF

#define ACCESS_ONCE(x) (*(volatile typeof(x) *)&(x))
#define READ_ONCE(x) ACCESS_ONCE(x)
#define WRITE_ONCE(x, val) ({ ACCESS_ONCE(x) = (val); })
#define MIN_INT(a, b) std::min<int>(a, b)
#define MAX_INT(a, b) std::max<int>(a, b)

#define MIN_LONG(a, b) std::min<long>(a, b)
#define MAX_LONG(a, b) std::max<long>(a, b)
#define LARGE_DRAM
#define MIN_DOUBLE(a, b) std::min<double>(a, b)
#define MAX_DOUBLE(a, b) std::max<double>(a, b)
#define LR_PRED(slope, minkey, key, fanout) std::min<int>(std::max<int>(static_cast<int>(slope * (key*1.0-minkey*1.0)), 0), fanout - 1)
//#define LR_PRED(b, slope, key, fanout) std::min<int>(std::max<int>(static_cast<int>(b + slope * key), 0), fanout - 1)
#define unlikely(x) __builtin_expect(!!(x),0)
#define likely(x) __builtin_expect(!!(x),1)
namespace plt {
    constexpr static uint32_t lockSet = ((uint32_t)1 << 31);
    constexpr static uint32_t lockMask = ((uint32_t)1 << 31) - 1;
    constexpr static uint8_t leafset = ((uint8_t)1 << 7) ;


    template<typename KeyType, typename ValueType>
    struct KV{
        KeyType first;
        ValueType second;
        bool operator==(const KV& other) const {
            return first == other.first ;
        }
    };
    struct __attribute__((aligned(64))) EntryFp {
        uint8_t fps[BUCKET_PER_ENTRY * BUCKET_LEN];
       // uintptr_t size:16;
        uintptr_t entryPtr;
        EntryFp *nextLevel;
        void clear() {
            memset(fps, 0, BUCKET_PER_ENTRY * BUCKET_LEN);
        }
    };
    struct BlockHeader {
        uint64_t ptr;
        uint8_t v;
        uint8_t num;
        uint8_t unused[2];
        uint32_t  version_lock;
        inline void get_lock(){
            uint32_t new_value = 0;
            uint32_t old_value = 0;
            do {
                while (true) {
                    old_value = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
                    if (!(old_value & lockSet)) {
                        break;
                    }
                }
                new_value = old_value | lockSet;
            } while (!CAS(&version_lock, &old_value, new_value));
        }
        inline bool try_get_lock() {
            uint32_t v = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
            if (v & lockSet) {
                return false;
            }
            auto old_value = v & lockMask;
            auto new_value = v | lockSet;
            return CAS(&version_lock, &old_value, new_value);
        }
        inline void release_lock() {
            uint32_t v = version_lock;
            __atomic_store_n(&version_lock, v + 1 - lockSet, __ATOMIC_RELEASE);
        }
        inline bool check_version_change(uint32_t old_version) const {
            auto value = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
            return (old_version != value);
        }
        inline bool check_size_change(uint32_t old_size) const {
            auto value = __atomic_load_n(&num, __ATOMIC_ACQUIRE);
            return (old_size != value);
        }

        inline void clearNum() {
            __atomic_store_n(&num, 0, __ATOMIC_RELEASE);
        }
        inline void addNum() {
            __atomic_store_n(&num, num + 1, __ATOMIC_RELEASE);
        }
        inline void release_lock_without_change() {
            __atomic_store_n(&version_lock, version_lock & lockMask, __ATOMIC_RELEASE);
        }
        void setBufferPtr(uint64_t bufptr){
            ptr=bufptr;
        }
        uint64_t getBufferPtr(){
            return ptr;
        }

    };
    struct __attribute__((aligned(32))) blockMeta{
    public:
        uint8_t v;
        uint8_t layerPtr[5];
        uint8_t layerfpPtr[6];
        uint32_t  version_lock;  // 1 bit lock & 31 bit block version
        uint8_t fps[15];
        uint8_t num;
        void clear() {
            memset(fps, 0, 15);
        }
        void setLayerPtr(uint64_t lptr){
           // layerPtr[0]=0;//static_cast<uint8_t>(lptr&0xff);
            layerPtr[0]=static_cast<uint8_t>(lptr>>8&0xff);
            layerPtr[1]=static_cast<uint8_t>((lptr>>16)&0xff);
            layerPtr[2]=static_cast<uint8_t>((lptr>>24)&0xff);
            layerPtr[3]=static_cast<uint8_t>((lptr>>32)&0xff);
            layerPtr[4]=static_cast<uint8_t>((lptr>>40)&0xff);
        }
        uint64_t getLayerPtr(){
            return //static_cast<uint64_t>(layerPtr[0])
                   static_cast<uint64_t>(layerPtr[0])<<8
                   |static_cast<uint64_t>(layerPtr[1])<<16
                   |static_cast<uint64_t>(layerPtr[2])<<24
                   |static_cast<uint64_t>(layerPtr[3])<<32
                   |static_cast<uint64_t>(layerPtr[4])<<40;
        }
        void setFPtr(uint64_t fptr){
            layerfpPtr[0]=static_cast<uint8_t>(fptr&0xff);
            layerfpPtr[1]=static_cast<uint8_t>(fptr>>8&0xff);
            layerfpPtr[2]=static_cast<uint8_t>((fptr>>16)&0xff);
            layerfpPtr[3]=static_cast<uint8_t>((fptr>>24)&0xff);
            layerfpPtr[4]=static_cast<uint8_t>((fptr>>32)&0xff);
            layerfpPtr[5]=static_cast<uint8_t>((fptr>>40)&0xff);
        }
        uint64_t getFPtr(){
            return static_cast<uint64_t>(layerfpPtr[0])
                   |static_cast<uint64_t>(layerfpPtr[1])<<8
                   |static_cast<uint64_t>(layerfpPtr[2])<<16
                   |static_cast<uint64_t>(layerfpPtr[3])<<24
                   |static_cast<uint64_t>(layerfpPtr[4])<<32
                   |static_cast<uint64_t>(layerfpPtr[5])<<40;
        }
        void freeFp(){
            uint64_t _fps=getFPtr();
            auto fp= reinterpret_cast<EntryFp*>(_fps);
            while(fp!= nullptr){
                auto next=fp[0].nextLevel;
#ifndef LARGE_DRAM
                free(fp);
#endif
                fp=next;
            }
        }
        inline void get_lock(){
            uint32_t new_value = 0;
            uint32_t old_value = 0;
            do {
                while (true) {
                    old_value = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
                    if (!(old_value & lockSet)) {
                        break;
                    }
                }
                new_value = old_value | lockSet;
            } while (!CAS(&version_lock, &old_value, new_value));
        }
        inline bool try_get_lock() {
            uint32_t v = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
            if (v & lockSet) {
                return false;
            }
            auto old_value = v & lockMask;
            auto new_value = v | lockSet;
            return CAS(&version_lock, &old_value, new_value);
        }
        inline void release_lock() {
            uint32_t v = version_lock;
            __atomic_store_n(&version_lock, v + 1 - lockSet, __ATOMIC_RELEASE);
        }

        inline void clearNum() {
            __atomic_store_n(&num, 0, __ATOMIC_RELEASE);
        }
        inline void addNum() {
            __atomic_store_n(&num, num + 1, __ATOMIC_RELEASE);
        }
        inline void release_lock_without_change() {
            __atomic_store_n(&version_lock, version_lock & lockMask, __ATOMIC_RELEASE);
        }
        // lock a block
        inline void get_lock(uint32_t global_version){
            //check_global_version(global_version);
            uint32_t new_value = 0;
            uint32_t old_value = 0;
            do {
                while (true) {
                    old_value = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
                    if (!(old_value & lockSet)) {
                        break;
                    }
                }
                new_value = old_value | lockSet;
            } while (!CAS(&version_lock, &old_value, new_value));
        }
        // check lock
        inline bool check_lock(uint32_t global_version) {
            // check_global_version(global_version);
            uint32_t v = __atomic_load_n(&version_lock,  __ATOMIC_ACQUIRE);
            if (v & lockSet) {
                return false;
            }
            return true;
        }
        // get the block version
        inline void get_version(uint32_t &version) const {
            do {
                version = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
            } while (version & lockSet);
        }
        // test whether the version has changed, if change, return true
        inline bool check_version_change(uint32_t old_version) const {
            auto value = __atomic_load_n(&version_lock, __ATOMIC_ACQUIRE);
            return (old_version != value);
        }
        inline bool check_size_change(uint32_t old_size) const {
            auto value = __atomic_load_n(&num, __ATOMIC_ACQUIRE);
            return (old_size != value);
        }
    };
    // A CDF coordinate.
    template<typename KeyType>
    struct Coord {
        KeyType x;
        double y;
    };

    template<typename KeyType>
    struct SegmentMessage {
        KeyType key;
        size_t offset;
        size_t size;
        float slope;
    };
    template<typename KeyType>
    std::vector<std::pair<KeyType*,int>>  treeLayout;
    template<typename KeyType>
    class innerSlot{
    public:
        float slope;
        uint32_t size;
        KeyType minkey;
        uint8_t v;
        uint8_t isLeaf;
        uint16_t childNodePtr[3];
        uint16_t bufferPtr[3];
        uint16_t levelnum; //
        innerSlot():minkey(0),slope(0),size(0),levelnum(0),v(0),isLeaf(0){}

        inline void addLevel() {
            __atomic_store_n(&levelnum, levelnum + 1, __ATOMIC_RELEASE);

        }
        void setSlot(KeyType _minkey,uintptr_t _ptr,float _slope,uint32_t _size,uint8_t _v,uint8_t _isLeaf){
            minkey=_minkey;
            setChildNodePtr(_ptr);
            slope=_slope;
            size=_size;
            v=_v;
            isLeaf=_isLeaf;
        }
        void setSlot(uintptr_t _ptr,uint8_t _v,uint8_t _isLeaf){
            setChildNodePtr(_ptr);
            v=_v;
            isLeaf=_isLeaf;
        }
        bool isleaf(){
            return (uint8_t)(isLeaf&leafset)==leafset;
        }
        void setBufferPtr(uint64_t ptr){
            bufferPtr[0]=static_cast<uint16_t>(ptr&0xffff);
            bufferPtr[1]=static_cast<uint16_t>((ptr>>16)&0xffff);
            bufferPtr[2]=static_cast<uint16_t>((ptr>>32)&0xffff);
        }
        void setChildNodePtr(uint64_t ptr){
            childNodePtr[0]=static_cast<uint16_t>(ptr&0xffff);
            childNodePtr[1]=static_cast<uint16_t>((ptr>>16)&0xffff);
            childNodePtr[2]=static_cast<uint16_t>((ptr>>32)&0xffff);
        }
        uint64_t getChildNodePtr(){
            return childNodePtr[0]|(static_cast<uint64_t>(childNodePtr[1])<<16)|(static_cast<uint64_t>(childNodePtr[2])<<32);

        }
        uint64_t getBufferPtr(){
            return bufferPtr[0]|(static_cast<uint64_t>(bufferPtr[1])<<16)|(static_cast<uint64_t>(bufferPtr[2])<<32);
        }
        inline void predict_block (const KeyType &key,int &predicted_block ,const int &kvsPerBlock)  {
            auto block=static_cast<float>(key-minkey) * slope /kvsPerBlock;
            predicted_block=block<0?0:(block>=size?size-1:static_cast<int>(block));
            return ;//predicted_block<0?0:(predicted_block>=size?size-1:static_cast<int>(predicted_block));
           // return std::max<int>(std::min<int>(predicted_block, size - 1), 0) ;
        }
        inline void check_global_version(uint8_t global_version){
            if (global_version != v) {
                __atomic_store_n(&isLeaf, leafset, __ATOMIC_RELEASE);
                __atomic_store_n(&v, global_version, __ATOMIC_RELEASE);
            }
        }
        inline bool try_get_split_lock(){
            uint8_t new_value = 0;
            uint8_t old_value = __atomic_load_n(&isLeaf, __ATOMIC_ACQUIRE);
            if (!(old_value & 1u)) {
                new_value = old_value | 1u;
                return CAS(&isLeaf, &old_value, new_value);
            }
            return false;
        }

        inline void get_write_lock(){
            ADD(&isLeaf,2);
        }

        inline void get_read_lock(){
            ADD(&isLeaf,4);
        }
        inline bool check_write_lock() {

            uint8_t v = __atomic_load_n(&isLeaf, __ATOMIC_ACQUIRE);
            if (v & 2u) {
                return false;
            }
            return true;
        }

        inline bool check_read_lock() {

            uint8_t v = __atomic_load_n(&isLeaf, __ATOMIC_ACQUIRE);
            if (v & 4u) {
                return false;
            }
            return true;
        }

        inline void release_lock() {
            __atomic_store_n(&isLeaf, leafset, __ATOMIC_RELEASE);
        }


    };
    template<typename KeyType>
    class innerNode{
    public:
        KeyType minkey;
        float slope;
        uint32_t fanout;
        uint8_t unused[16];
        innerSlot<KeyType> slots[0];
        innerNode():minkey(0),fanout(0),slope(0){}
        void setInnerNode(KeyType* keys,uint32_t n_keys){
            assert(n_keys >= 0);
            minkey = keys[0];
            if (n_keys >= 2) {
                slope= 1.0 * n_keys / (keys[n_keys - 1] - keys[0]);
                fanout = n_keys ;
            } else {
                slope = 0;
                fanout = 1;
            }
        }
        void setInnerNode(KeyType min,KeyType max,int _size){
            minkey=min;
            fanout=_size;
            slope=_size*1.0/(max-min);

        }
/*        innerSlot<KeyType> * findSlot(KeyType key){
            int idx=LR_PRED(b,slope,key,fanout);
           return &slots[idx];
        }*/
    };

    template<typename KeyType, typename ValueType>
    struct __attribute__((aligned(256))) Block{
        typedef  struct KV<KeyType,ValueType> KV;
        BlockHeader head;
        KV kv[(BLOCK_SIZE-sizeof(BlockHeader))/sizeof(KV)];
    };

    template<typename KeyType>
    inline int predict_block (KeyType key, innerSlot<KeyType> *  seginf)  {
        int predicted_block = (key-seginf->minkey) * seginf->slope /15;
        return std::max<int>(std::min<int>(predicted_block, seginf->size - 1), 0) ;
    }
    template<typename KeyType>
    inline int predict_block (const KeyType &key, const KeyType &minkey,const float &slope ,const uint32_t &size,const int &kvsPerBlock)  {
        int predicted_block = (key-minkey) * slope /kvsPerBlock;
        return std::max<int>(std::min<int>(predicted_block, size - 1), 0) ;
    }
    template<typename KeyType>
    inline int array_lower_bound(KeyType* data, double x, int l, int r) {
        while (l < r) {
            int mid = (l + r) >> 1;
            if (data[mid] >= x) {
                r = mid;
            } else {
                l = mid + 1;
            }
        }
        return l;
    }


    template<typename KeyType, typename ValueType>
    inline void linearReg_w_simple_strategy(const std::vector<KV<KeyType,ValueType>> &X, int start, double &a, double &b, uint32_t n) {
        int left_n = n / 2;
        int right_n = n - left_n - 1;
        //left_n+=start;
        //right_n+=start;
        KeyType x_middle = X[left_n+start].first;
        double left_slope = 1.0 * left_n / (x_middle*1.0 - X[start].first*1.0);
        double right_slope = 1.0 * right_n / (X[start+n-1].first*1.0 - x_middle*1.0);
        b = MAX_DOUBLE(left_slope, right_slope);
        a = -b * x_middle + left_n;
    }

    template<typename KeyType, typename ValueType>
    inline void linearReg_(const std::vector<KV<KeyType,ValueType>> &X, int start,double &a, double &b, uint32_t n) {
        double nu_b = 0;
        double de_b = 0;
        double mean_x = 0;
        int end= start + n;
        for (int i = start; i < end; ++i) {
            de_b += X[i].first * 1.0 * i;
            nu_b += X[i].first * 1.0 * X[i].first;
            mean_x += X[i].first;
        }

        mean_x /= n;
        double mean_y = (n - 1) / 2.0;
        de_b -= mean_x * mean_y * n;
        nu_b -= mean_x * mean_x * n;

        if (nu_b != 0) {
            b = de_b / nu_b;
        }  else {
            b = 0;
        }
        a = mean_y - b * mean_x;
    }

    template<typename KeyType, typename ValueType>
    inline void linearReg_w_expanding(const std::vector<KV<KeyType,ValueType>> & X,int start, double &a, double &b,  uint32_t n, bool use_simple_strategy) {
        if (!use_simple_strategy) {
            linearReg_(X, start,a, b, n);
        } else {
            linearReg_w_simple_strategy(X,start, a, b, n);
        }

    }

    struct SearchBound {
        size_t begin;
        size_t end; // Exclusive.
    };



    /*** hash function ***/
    template <class T>
    static inline unsigned char hashcode1B(T y) {
        uint64_t x = (uint64_t &)y;
        x ^= x>>32;
        x ^= x>>16;
        x ^= x>>8;
        return (unsigned char)(x&0x0ffULL);
    }




} // namespace plt

#endif //ART_TEST_COMMON_H
