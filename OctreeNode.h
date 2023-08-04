#include <vector>
#ifndef POINT_H
#define POINT_H
// Octree node structure

template <typename T>
struct OctreeNodeBase : T
{
    bool is_leaf = true;
    void* children = nullptr;
};

template <typename T>
struct OctreeLeafNode{
    T value;
};

template<typename T>
using OctreeInnerNode = OctreeNodeBase<OctreeLeafNode<T>>;

#endif