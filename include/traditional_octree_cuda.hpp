// Define a basic point structure
#include <iostream>
#include <vector>
#include "cuda/cuda_device.h"
#include "traditional_octree.hpp"

// Octree class
class OctreeCuda {
public:
    redwood::lang::cuda::CudaDevice* _device;  
    OctreeNode* _root;
    std::vector<redwood::lang::DeviceAllocation> buffers;

    OctreeCuda(const BoundingBox& bounds, redwood::lang::cuda::CudaDevice* device) {
        _device = device;
        redwood::lang::DeviceAllocation root_buf;
        _device->allocate_memory(redwood::lang::Device::AllocParams{sizeof(OctreeNode)}, &root_buf);
        void* mappedData;
        buffers.push_back(root_buf);
        root_buf.device->map(root_buf, &mappedData);
        _root = static_cast<OctreeNode*>(mappedData);
    }

    ~OctreeCuda() {
        for(int i = 0; i < buffers.size(); ++i){
            _device->dealloc_memory(buffers[i]);
        }
    }

    // Insert a point into the octree
    void insert(const Point& point) {
        if (_root->bounds.contains(point)) {
            insert(_root, point);
        }
    }

    // Recursive insert method
    void insert(OctreeNode* node, const Point& point) {
        if (node->isLeaf()) {
            node->points.push_back(point);
            if (node->points.size() > 8) {
                subdivide(node);
            }
        } else {
            redwood::lang::DeviceAllocation child_buff;
            _device->allocate_memory(redwood::lang::Device::AllocParams{sizeof(OctreeNode)}, &child_buff);
            void* mappedData;
            child_buff.device->map(child_buff, &mappedData);
            buffers.push_back(child_buff);
            OctreeNode* node = static_cast<OctreeNode*>(mappedData);
            int index = getSubnodeIndex(node, point);
            if (node->children[index] == nullptr) {
                BoundingBox childBounds = getChildBounds(node->bounds, index);
                *node = OctreeNode(childBounds);
                node->children[index] = node;
            }
            insert(node->children[index], point);
        }
    }

    // Subdivide the node into 8 subnodes
    void subdivide(OctreeNode* node) {
        const BoundingBox& currentBounds = node->bounds;
        redwood::lang::DeviceAllocation child_buff;
         _device->allocate_memory(redwood::lang::Device::AllocParams{sizeof(OctreeNode)*8}, &child_buff);
        void* mappedData;
        child_buff.device->map(child_buff, &mappedData);
        buffers.push_back(child_buff);
        OctreeNode* childrens = static_cast<OctreeNode*>(mappedData);
        for (int i = 0; i < 8; ++i) {
            BoundingBox childBounds = getChildBounds(currentBounds, i);
            childrens[i] = OctreeNode(childBounds);
            node->children[i] = &childrens[i];
        }

        // Move the points to the subnodes
        for (const Point& point : node->points) {
            int index = getSubnodeIndex(node, point);
            insert(node->children[index], point);
        }
        node->points.clear();
    }

    // Get the index of the subnode for a point
    int getSubnodeIndex(OctreeNode* node, const Point& point) {
        const BoundingBox& currentBounds = node->bounds;
        Point midpoint((currentBounds.min.x + currentBounds.max.x) / 2,
                       (currentBounds.min.y + currentBounds.max.y) / 2,
                       (currentBounds.min.z + currentBounds.max.z) / 2);

        int index = 0;
        if (point.x > midpoint.x) index |= 1;
        if (point.y > midpoint.y) index |= 2;
        if (point.z > midpoint.z) index |= 4;
        return index;
    }

    // Calculate the bounding box for a subnode
    BoundingBox getChildBounds(const BoundingBox& currentBounds, int index) {
        Point min = currentBounds.min;
        Point max = currentBounds.max;
        Point midpoint((min.x + max.x) / 2, (min.y + max.y) / 2, (min.z + max.z) / 2);

        if (index & 1) min.x = midpoint.x;
        else max.x = midpoint.x;
        if (index & 2) min.y = midpoint.y;
        else max.y = midpoint.y;
        if (index & 4) min.z = midpoint.z;
        else max.z = midpoint.z;

        return BoundingBox(min, max);
    }
};
