#include <vector>
#include "Point.hpp"
// Octree node structure
struct OctreeNode {
    std::vector<Point> points;
    int children[8]; // Array to store the child indices

    OctreeNode() {
        for (int i = 0; i < 8; ++i) {
            children[i] = -1; // Initialize children to -1 (invalid index)
        }
    }
};