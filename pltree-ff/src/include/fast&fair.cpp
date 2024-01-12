//
// Created by mcooiedo on 23-9-21.
//
#include "fast&fair.h"

bool btree::find(_key_t key, _payload_t & val){
    page* p = root;

    while(p->hdr.leftmost_ptr != NULL) {
        p = (page *)(p->linear_search(key));
    }

    page *t;
    while((t = (page *)p->linear_search(key)) == p->hdr.sibling_ptr) {
        p = t;
        if(!p) {
            break;
        }
    }

    if(!t) {
        //printf("NOT FOUND %lu, t = %x\n", key, t);
        return false;
    }

    val = (_payload_t)t;
    return true;
}

// insert the key in the leaf node
void btree::insert(_key_t key, _payload_t right){ // need to be std::string
    page* p = root;

    while(p->hdr.leftmost_ptr != NULL) {
        p = (page *)(p->linear_search(key));
    }
    if(!p->store(this, NULL, key, (char *)right, true)) { // store
        insert(key, right);
    }
}

// store the key into the node at the given level
void btree::btree_insert_internal(char *left, _key_t key, char *right, uint32_t level) {
    page* p = root;

    if(level > p->hdr.level)
        return;

    while(p->hdr.level > level)
        p = (page *)(p->linear_search(key));

    if(!p->store(this, NULL, key, right, true)) {
        btree_insert_internal(left, key, right, level);
    }
}

bool btree::remove(_key_t key) {
    page* p = root;

    while(p->hdr.leftmost_ptr != NULL){
        p = (page *)(p->linear_search(key));
    }

    page *t;
    while((t = (page *)p->linear_search(key)) == p->hdr.sibling_ptr) {
        p = t;
        if(!p)
            break;
    }

    if(p) {
        //Modified by ypluo, eliminate the dead loop on delete key
        return p->remove_rebalancing(this, key, false);

        if(p->hdr.is_deleted == 1) {
            //galc->free(p);
        }

    } else {
        //printf("not found the key to delete %lu\n", key);
        return false;
    }
}

void btree::btree_delete_internal(_key_t key, char *ptr, uint32_t level, _key_t *deleted_key,
                                  bool *is_leftmost_node, page **left_sibling) {
    page* p = root;

    if(level > p->hdr.level)
        return;


    while(p->hdr.level > level) {
        p = (page *)(p->linear_search(key));
    }



    if((char *)p->hdr.leftmost_ptr == (char *)(ptr)) {
        *is_leftmost_node = true;
        return;
    }

    *is_leftmost_node = false;

    for(int i=0; p->records[i].val != NULL; ++i) {
        if(p->records[i].val == (char *)(ptr)) {
            if(i == 0) {
                if((char *)p->hdr.leftmost_ptr != p->records[i].val) {
                    *deleted_key = p->records[i].key;
                    *left_sibling = (page *)(p->hdr.leftmost_ptr);

                    p->remove_rebalancing(this, *deleted_key, false); // Modified by ypluo
                    //p->remove(this, *deleted_key, false, false);
                    break;
                }
            }
            else {
                if(p->records[i - 1].val != p->records[i].val) {
                    *deleted_key = p->records[i].key;
                    *left_sibling = (page *)(p->records[i - 1].val);
                    p->remove_rebalancing(this, *deleted_key, false); // Modified by ypluo
                    //p->remove(this, *deleted_key, false, false);
                    break;
                }
            }
        }
    }
}

// insert the key in the leaf node
bool btree::update(_key_t key, _payload_t value){ //need to be std::string
    page* p = root;

    while(p->hdr.leftmost_ptr != NULL) {
        p = (page *)(p->linear_search(key));
    }

    return p->update_value(key, value);
}

// Upsert
// 2 : Goto learned node to insert; 3 : update; 4 : insert
uint32_t btree::upsert(_key_t key, _payload_t value, uint32_t ds){
    page* p = root;

    while(p->hdr.leftmost_ptr != NULL) {
        p = (page *)(p->linear_search(key));
    }
    page* first_p = p;
    page *t;
    while((t = (page *)p->linear_search(key)) == p->hdr.sibling_ptr) {
        p = t;
        if(!p) {
            break;
        }
    }

    if(!t) {
        if (ds > 0) {
            return 2;
        }
        if(!first_p->store(this, NULL, key, (char *)value, true)) {
            insert(key, value);
        }
        return 4;
    }

    first_p->update_value(key, value);
    return 3;

}

// Function to search keys from "min" to "max"
int btree::scan(_key_t min, int scan_sz, _payload_t * buf) {
    page* p = root;

    while(p) {
        if(p->hdr.leftmost_ptr != NULL) {
            // The current page is internal
            p = (page *)(p->linear_search(min));
            //p = (page *)p->linear_search(min);
        }
        else {
            // Found a leaf
            return p->linear_search_range(min, scan_sz, buf);
        }
    }

    return 0;
}

void btree::range_query (_key_t lower_bound, _key_t upper_bound, std::vector<std::pair<_key_t, _payload_t>>& answers) {
    page* p = root;

    while (p->hdr.leftmost_ptr != NULL) {
        p = (page *)(p->hdr.leftmost_ptr);
    }

    int flag = p->range_query(lower_bound, upper_bound, answers);
    while (flag == 1 && p->hdr.sibling_ptr) {
        p = (page *)(p->hdr.sibling_ptr);
        flag = p->range_query(lower_bound, upper_bound, answers);
    }
}
/*
int btree::range_query_by_size (_key_t lower_bound, int toscan, std::map<_key_t, _payload_t>& answers) {
    page* p = root;

    while (p->hdr.leftmost_ptr != NULL) {
        p = (page *)(p->hdr.leftmost_ptr);
    }

    int num = p->range_query(lower_bound, toscan, answers);
    while (num >= 1 && p->hdr.sibling_ptr) {
        p = (page *)(p->hdr.sibling_ptr);
        num = p->range_query(lower_bound, num, answers);
    }
    return num;
}
*/


void btree::get_data (std::vector<_key_t>& keys, std::vector<_payload_t>& payloads) {
    page* p = root;

    while (p->hdr.leftmost_ptr != NULL) {
        p = (page *)(p->hdr.leftmost_ptr);
    }

    p->get_data(keys, payloads);
    while (p->hdr.sibling_ptr) {
        p = (page *)(p->hdr.sibling_ptr);
        p->get_data(keys, payloads);
    }

}
/*void btree::destroy(){
    page *leftmost = root;
    do {
        page *sibling = leftmost;
        while(sibling) {
            if(sibling->hdr.level == 0) {
                sibling->destroy(this);
            }
            sibling = (page *) (sibling->hdr.sibling_ptr);
        }
        leftmost = (page *) (leftmost->hdr.leftmost_ptr);
    } while(leftmost);
}*/

void btree::node_size (uint64_t& overflow_number) {

    page *leftmost = root;
    do {
        page *sibling = leftmost;
        while(sibling) {
            if(sibling->hdr.level == 0) {
                overflow_number += sibling->hdr.last_index + 1;
            }
            sibling = (page *) (sibling->hdr.sibling_ptr);
        }
        leftmost = (page *) (leftmost->hdr.leftmost_ptr);
    } while(leftmost);

}

void btree::printAll(){
    int total_keys = 0;
    page *leftmost = root;

    printf("root: %lx\n", (long unsigned int)leftmost);
    do {
        page *sibling = leftmost;
        while(sibling) {
            if(sibling->hdr.level == 0) {
                total_keys += sibling->hdr.last_index + 1;
            }
            sibling->print();
            sibling = (page *) (sibling->hdr.sibling_ptr);
        }
        printf("-----------------------------------------\n");
        leftmost = (page *) (leftmost->hdr.leftmost_ptr);
    } while(leftmost);

    printf("total number of keys: %d\n", total_keys);
}

// additional function to iterater the tree
char ** btree::find_lower(_key_t key) const{
    page* p = root;

    while(p->hdr.leftmost_ptr != NULL) {
        p = (page *)(p->linear_search(key));
    }

    page *t = NULL;
    _key_t discard_key;
    while((t = (page *)p->linear_search_lower(key, discard_key)) == p->hdr.sibling_ptr) {
        p = (page *)(t);
        if(!p) {
            break;
        }
    }
    return (char **)t;
}

bool btree::try_remove(_key_t key, bool & need_rebuild) {
    page* p = root;
    while(p->hdr.leftmost_ptr != NULL){
        p = (page *)(p->linear_search(key));
    }
    page *t;
    _key_t new_key = -1;
    while((t = (page *)p->linear_search_lower(key, new_key)) == p->hdr.sibling_ptr) {
        p = (page *)(t);
        if(!p)
            break;
    }

    p->remove_rebalancing(this, new_key);
    if(p->hdr.is_deleted == 1) {
        //galc->free(p);
    }

    need_rebuild = false;
    return true;
    return new_key;
}

void btree::print_height() {
    page* p = root;
    printf("Height: %d\n", p->hdr.level + 1);
}

// btree::btree() {
//     root = (new page());
// }

btree::btree(void** addr, bool init)
        : addr(addr) {
    if (init) {
        root = new page();
        *addr =(root);
        Alloc::do_flush_with_double_fence(addr, 8);
    }
    else {
        root = (page*) (*addr);
    }
}

void btree::setNewRoot(page *new_root) {
    root = new_root;
    *addr =(new_root);
    Alloc::do_flush_with_double_fence(addr, 8);
}
