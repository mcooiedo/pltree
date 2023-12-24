//
// Created by Cshuang on 2020/10/12.
//

#include "art_tree.h"
#include "segment.h"

template<> uint64_t palert::ArtTree<uint32_t>::allocated_byte_count = 0;
template<> uint64_t palert::ArtTree<uint64_t>::allocated_byte_count = 0;

template<> uint64_t palert::ArtTree<uint64_t>::node4_num = 0;
template<> uint64_t palert::ArtTree<uint64_t>::node16_num = 0;
template<> uint64_t palert::ArtTree<uint64_t>::node48_num = 0;
template<> uint64_t palert::ArtTree<uint64_t>::node256_num = 0;

/*template<> uint64_t  palert::Segment<uint32_t, uint32_t>::segment_allocated_byte = 0;
template<> uint64_t  palert::Segment<uint64_t, uint64_t>::segment_allocated_byte = 0;*/
