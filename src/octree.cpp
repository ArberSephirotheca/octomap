#include "octree.hpp"

#include <iostream>
#include <thread>

#include "binary_radix_tree.hpp"
namespace oct {

template<typename DATA_TYPE>
void OctNode<DATA_TYPE>::SetChild(const int child, const int my_child_idx) {
  children[my_child_idx] = child;
  // TODO: atomicOr in CUDA
  child_node_mask |= (1 << my_child_idx);
}

template<typename DATA_TYPE>
void OctNode<DATA_TYPE>::SetLeaf(const int leaf, const int my_child_idx) {
  children[my_child_idx] = leaf;
  // TODO: atomicOr in CUDA
  child_leaf_mask |= (1 << my_child_idx);
}

template<typename DATA_TYPE>
  OctNode<DATA_TYPE>* Octree<DATA_TYPE>::search(const Code_t prefix, const int code_len, const int oct_idx) {
  const OctNode<DATA_TYPE>& node = nodes_[oct_idx];
  for (int i = 0; i < 8; ++i) {
    Code_t new_pref = (prefix << 3) | i;
    if (node.child_node_mask & (1 << i)) {
      search(new_pref, code_len + 3, node.children[i]);
    }
    if (node.child_leaf_mask & (1 << i)) {
      return node.children[i];
    }
  }
}

}  // namespace oct