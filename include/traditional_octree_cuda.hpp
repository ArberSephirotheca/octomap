// Define a basic point structure
#pragma once
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



// Define a basic point structure
#pragma once
#include <iostream>
#include <vector>
#include "cpu/cpu_device.h"
#include "traditional_octree_base.hpp"

class OctreeNodeCUDA : public OctreeNodeBase
{

public:
	Point point;
	redwood::lang::DeviceAllocation *handle;

	// Represent the boundary of the cube
	OctreeNodeCUDA *children[8];
	Point topLeftFront, bottomRightBack;
	OctreeNodeCUDA()
	{
		point = Point();
	}
	OctreeNodeCUDA(int x, int y, int z)
	{
		point = Point(x, y, z);
	}
	OctreeNodeCUDA(int x1, int y1, int z1, int x2, int y2, int z2)
	{	
		// TODO: This is a bug

		// This use to construct Octree
		// with boundaries defined
		if (x2 < x1 || y2 < y1 || z2 < z1)
		{
			return;
		}
		point.x = -2;
		topLeftFront = Point(x1, y1, z1);
		bottomRightBack = Point(x2, y2, z2);
		for (int i = 0; i < 8; ++i)
		{
			children[i] = nullptr;
		}
		for (int i = TopLeftFront;
			 i <= BottomLeftBack;
			 ++i)
			*children[i] = OctreeNodeCUDA();
	}
	bool find(int x, int y, int z)
	{
		// If point is out of bound
		if (x < topLeftFront.x || x > bottomRightBack.x || y < topLeftFront.y || y > bottomRightBack.y || z < topLeftFront.z || z > bottomRightBack.z)
			return 0;

		// Otherwise perform binary search
		// for each ordinate
		int midx = (topLeftFront.x + bottomRightBack.x) / 2;
		int midy = (topLeftFront.y + bottomRightBack.y) / 2;
		int midz = (topLeftFront.z + bottomRightBack.z) / 2;

		int pos = -1;

		// Deciding the position
		// where to move
		if (x <= midx)
		{
			if (y <= midy)
			{
				if (z <= midz)
					pos = TopLeftFront;
				else
					pos = TopLeftBottom;
			}
			else
			{
				if (z <= midz)
					pos = BottomLeftFront;
				else
					pos = BottomLeftBack;
			}
		}
		else
		{
			if (y <= midy)
			{
				if (z <= midz)
					pos = TopRightFront;
				else
					pos = TopRightBottom;
			}
			else
			{
				if (z <= midz)
					pos = BottomRightFront;
				else
					pos = BottomRightBack;
			}
		}

		// If an internal node is encountered
		if (children[pos]->point.x == -2)
		{
			return children[pos]->find(x, y, z);
		}
		// If an empty node is encountered
		if (children[pos]->point.x == -1)
		{
			return false;
		}
		else
		{
			// If node is found with
			// the given value
			if (x == children[pos]->point.x && y == children[pos]->point.y && z == children[pos]->point.z)
				return 1;
		}
		return 0;
	}
};

// Octree class
class OctreeCUDA : public OctreeBase
{
public:
	redwood::lang::cuda::CudaDevice* _device;
	OctreeNodeCUDA *_root;
	std::vector<redwood::lang::DeviceAllocation *> buffers;

	OctreeCUDA(int x1, int y1, int z1, int x2, int y2, int z2, redwood::lang::cpu::CpuDevice *device)
	{
		// TODO: This is a bug
		std::cout << "allocate" << std::endl;
		_device = device;
		redwood::lang::DeviceAllocation *root_buf = new redwood::lang::DeviceAllocation();
		_device->allocate_memory(redwood::lang::Device::AllocParams{sizeof(OctreeNode)}, &root_buf);
		void *mappedData;
		buffers.push_back(root_buf);
		root_buf->device->map(*root_buf, &mappedData);
		_root = static_cast<OctreeNodeCUDA *>(mappedData);
        _device->upload_data(*root_buf, _root, sizeof(OctreeNodeCUDA, 1));
		* _root = OctreeNodeCUDA(x1, y1, z1);
		_root->handle = root_buf;

		std::cout<<"done allocate"<<std::endl;
	}

	~OctreeCUDA()
	{
		for (int i = 0; i < buffers.size(); ++i)
		{
			_device->dealloc_memory(*buffers[i]);
		}
	}

	void insert(int x, int y, int z)
	{
		// If the point already exists in the octree
		if (find(x, y, z) == true)
		{
			return;
		}
		insert(_root, x, y, z);
	}
	bool find(int x, int y, int z)
	{
		return _root->find(x, y, z);
	}
	void insert(OctreeNodeCUDA *node, int x, int y, int z)
	{
		// If the point is out of bounds
		if (x < node->topLeftFront.x || x > node->bottomRightBack.x || y < node->topLeftFront.y || y > node->bottomRightBack.y || z < node->topLeftFront.z || z > node->bottomRightBack.z)
		{
			return;
		}

		// Binary search to insert the point
		int midx = (node->topLeftFront.x + node->bottomRightBack.x) / 2;
		int midy = (node->topLeftFront.y + node->bottomRightBack.y) / 2;
		int midz = (node->topLeftFront.z + node->bottomRightBack.z) / 2;

		int pos = -1;

		// Checking the octant of
		// the point
		if (x <= midx)
		{
			if (y <= midy)
			{
				if (z <= midz)
					pos = TopLeftFront;
				else
					pos = TopLeftBottom;
			}
			else
			{
				if (z <= midz)
					pos = BottomLeftFront;
				else
					pos = BottomLeftBack;
			}
		}
		else
		{
			if (y <= midy)
			{
				if (z <= midz)
					pos = TopRightFront;
				else
					pos = TopRightBottom;
			}
			else
			{
				if (z <= midz)
					pos = BottomRightFront;
				else
					pos = BottomRightBack;
			}
		}

		// If an internal node is encountered
		if (node->children[pos]->point.x == -2)
		{
			 std::cout<<"internal point"<<std::endl;
			insert(node->children[pos], x, y, z);
			return;
		}

		// If an empty node is encountered
		else if (node->children[pos]->point.x == -1)
		{
			 std::cout <<"empty point"<<std::endl;
			_device->dealloc_memory(*node->children[pos]->handle);
			redwood::lang::DeviceAllocation *child_buf = new redwood::lang::DeviceAllocation();
			_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCUDA)), true, true}, child_buf);
			void *mappedData;
			buffers.push_back(child_buf);
			child_buf->device->map(*child_buf, &mappedData);
			OctreeNodeCUDA *node = static_cast<OctreeNodeCUDA *>(mappedData);
			*node = OctreeNodeCUDA(x, y, z);
			node->handle = child_buf;
			node->children[pos] = node;
			return;
		}
		else
		{
			std::cout<<"point"<<std::endl;
			int x_ = node->children[pos]->point.x,
				  y_ = node->children[pos]->point.y,
				  z_ = node->children[pos]->point.z;
			_device->dealloc_memory(*node->children[pos]->handle);
			node->children[pos] = nullptr;
			if (pos == TopLeftFront)
			{
				redwood::lang::DeviceAllocation *child_buf = new redwood::lang::DeviceAllocation();
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCUDA)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				child_buf->device->map(*child_buf, &mappedData);
				OctreeNodeCUDA *node = static_cast<OctreeNodeCUDA *>(mappedData);
				*node = OctreeNodeCUDA(node->topLeftFront.x,
									  node->topLeftFront.y,
									  node->topLeftFront.z,
									  midx,
									  midy,
									  midz);
				node->children[pos] = node;
			}

			else if (pos == TopRightFront)
			{
				redwood::lang::DeviceAllocation *child_buf = new redwood::lang::DeviceAllocation();
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCUDA)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				child_buf->device->map(*child_buf, &mappedData);
				OctreeNodeCUDA *node = static_cast<OctreeNodeCUDA *>(mappedData);
				*node = OctreeNodeCUDA(midx + 1,
									  node->topLeftFront.y,
									  node->topLeftFront.z,
									  node->bottomRightBack.x,
									  midy,
									  midz);
				node->children[pos] = node;
			}
			else if (pos == BottomRightFront)
			{
				redwood::lang::DeviceAllocation *child_buf = new redwood::lang::DeviceAllocation();
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCUDA)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				child_buf->device->map(*child_buf, &mappedData);
				OctreeNodeCUDA *node = static_cast<OctreeNodeCUDA *>(mappedData);
				*node = OctreeNodeCUDA(midx + 1,
									  midy + 1,
									  node->topLeftFront.z,
									  node->bottomRightBack.x,
									  node->bottomRightBack.y,
									  midz);
				node->children[pos] = node;
			}
			else if (pos == BottomLeftFront)
			{
				redwood::lang::DeviceAllocation *child_buf = new redwood::lang::DeviceAllocation();
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCUDA)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				child_buf->device->map(*child_buf, &mappedData);
				OctreeNodeCUDA *node = static_cast<OctreeNodeCUDA *>(mappedData);
				*node = OctreeNodeCUDA(node->topLeftFront.x,
									  midy + 1,
									  node->topLeftFront.z,
									  midx,
									  node->bottomRightBack.y,
									  midz);
				node->children[pos] = node;
			}
			else if (pos == TopLeftBottom)
			{
				redwood::lang::DeviceAllocation *child_buf = new redwood::lang::DeviceAllocation();
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCUDA)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				child_buf->device->map(*child_buf, &mappedData);
				OctreeNodeCUDA *node = static_cast<OctreeNodeCUDA *>(mappedData);
				*node = OctreeNodeCUDA(node->topLeftFront.x,
									  node->topLeftFront.y,
									  midz + 1,
									  midx,
									  midy,
									  node->bottomRightBack.z);
				node->children[pos] = node;
			}
			else if (pos == TopRightBottom)
			{
				redwood::lang::DeviceAllocation *child_buf = new redwood::lang::DeviceAllocation();
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCUDA)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				child_buf->device->map(*child_buf, &mappedData);
				OctreeNodeCUDA *node = static_cast<OctreeNodeCUDA *>(mappedData);
				*node = OctreeNodeCUDA(midx + 1,
									  node->topLeftFront.y,
									  midz + 1,
									  node->bottomRightBack.x,
									  midy,
									  node->bottomRightBack.z);
				node->children[pos] = node;
			}
			else if (pos == BottomRightBack)
			{
				redwood::lang::DeviceAllocation *child_buf = new redwood::lang::DeviceAllocation();
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCUDA)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				child_buf->device->map(*child_buf, &mappedData);
				OctreeNodeCUDA *node = static_cast<OctreeNodeCUDA *>(mappedData);
				*node = OctreeNodeCUDA(midx + 1,
									  midy + 1,
									  midz + 1,
									  node->bottomRightBack.x,
									  node->bottomRightBack.y,
									  node->bottomRightBack.z);
				node->children[pos] = node;
			}
			else if (pos == BottomLeftBack)
			{
				redwood::lang::DeviceAllocation *child_buf = new redwood::lang::DeviceAllocation();
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCUDA)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				child_buf->device->map(*child_buf, &mappedData);
				OctreeNodeCUDA *node = static_cast<OctreeNodeCUDA *>(mappedData);
				*node = OctreeNodeCUDA(node->topLeftFront.x,
									  midy + 1,
									  midz + 1,
									  midx,
									  node->bottomRightBack.y,
									  node->bottomRightBack.z);
				node->children[pos] = node;
			}
			insert(node->children[pos], x_, y_, z_);
			insert(node->children[pos], x, y, z);
		}
	}
};