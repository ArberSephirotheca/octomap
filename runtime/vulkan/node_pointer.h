#pragma once
#include <atomic>
#include "struct/snode.h"

namespace redwood::lang{
class PointerNode: public SNode {
public:
PointerNode(int depth,
        SNodeType t)
    : SNode(depth, t){
    _ = false;
    for(int i = 0; i < redwood_max_num_elements; i++){
        locks[i] = 0;
        data_pointers[i] = nullptr;
    }
    RW_INFO("Created PointerNode");
}

bool is_representative(uint32_t mask, uint64_t value);

 int Pointer_get_num_elements() {
  return max_num_elements();
}
void Pointer_activate(int i);

void Pointer_deactivate(int i);

bool Pointer_is_active(int i);

Ptr Pointer_lookup_element(int i);

//void set_cell_size_bytes(std::size_t size) override;
private:
    bool _ = false;
    std::atomic<uint64_t> locks[redwood_max_num_elements];
    Ptr data_pointers[redwood_max_num_elements];
};
}  // namespace rewood::lang


