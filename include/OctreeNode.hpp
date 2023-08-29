#include <vector>
#ifndef POINT_H
#define POINT_H
// Octree node structure

template <typename T>
struct OctreeNodeBase : T
{
    bool is_leaf = true;
    int* children = nullptr;
    uint8_t child_node_mask;
    uint8_t child_leaf_mask;
};

template <typename T>
struct OctreeLeafNode{
    T value;
};

template<typename T>
using OctreeInnerNode = OctreeNodeBase<OctreeLeafNode<T>>;

#endif