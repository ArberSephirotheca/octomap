#pragma once
#include "struct/snode.h"

namespace redwood::lang{
class DynamicNode: public SNode {
public:
DynamicNode(int depth,
        SNodeType t)
    : SNode(depth, t){
  ptr = nullptr;
  lock = false;
  n = 0;
}
void Dynamic_activate(int i);


void Dynamic_deactivate();

Ptr Dynamic_allocate(int32_t *len);

bool Dynamic_is_active(int i) {
  return i < n;
}

Ptr Dynamic_lookup_element( int i);

int32_t Dynamic_get_num_elements() {
  return n;
}

//void set_cell_size_bytes(std::size_t size) override;
private:
  std::atomic<bool> lock;
  // num of elements
  std::atomic<int> n;
  // pointer to chunk
  Ptr ptr;
};
}  // namespace rewood::lang