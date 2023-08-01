#include <vector>
#ifndef POINT_HPP
#define POINT_HPP
// Octree node structure
struct OctreeNode {
    uint32_t code;  // Morton code for the node's position
    int depth;      // Depth level in the octree

    int data_index; // Index to the first data point in the node
    int children[8]; // Array to store the child indices

    OctreeNode(uint32_t code, int depth) {
        for (int i = 0; i < 8; ++i) {
            children[i] = -1; // Initialize children to -1 (invalid index)
        }
    }
};

#endif