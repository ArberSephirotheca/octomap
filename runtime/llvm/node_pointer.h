#pragma once
#include "struct/snode.h"

namespace redwood::lang{
class PointerNode: public SNode {
public:
PointerNode(int depth,
        SNodeType t)
    : SNode(depth, t){
    _ = false;
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
};
}  // namespace rewood::lang


