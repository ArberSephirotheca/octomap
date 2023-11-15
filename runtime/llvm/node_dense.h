#pragma once

// Specialized Attributes and functions
struct DenseMeta : public StructMeta {
  int morton_dim;
};

using Ptr = uint8_t*;
int32_t Dense_get_num_elements(Ptr meta, Ptr node) {
  return ((StructMeta *)meta)->max_num_elements;
}

void Dense_activate(Ptr meta, Ptr node, int i) {
  // Dense elements are always active
}

bool Dense_is_active(Ptr meta, Ptr node, int i) {
  return true;
}

Ptr Dense_lookup_element(Ptr meta, Ptr node, int i) {
  return node + ((StructMeta *)meta)->element_size * i;
}