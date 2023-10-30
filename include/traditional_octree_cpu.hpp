// Define a basic point structure
#pragma once
#include <iostream>
#include <vector>
#include "cpu/cpu_device.h"
#include "traditional_octree_base.hpp"

class OctreeNodeCPU : public OctreeNodeBase
{

public:
	Point point;
	redwood::lang::DeviceAllocation *handle;

	// Represent the boundary of the cube
	OctreeNodeCPU *children[8];
	Point topLeftFront, bottomRightBack;
	OctreeNodeCPU()
	{
		point = Point();
	}
	OctreeNodeCPU(float x, float y, float z)
	{
		point = Point(x, y, z);
	}
	OctreeNodeCPU(float x1, float y1, float z1, float x2, float y2, float z2)
	{
		// This use to construct Octree
		// with boundaries defined
		if (x2 < x1 || y2 < y1 || z2 < z1)
		{
			return;
		}
		point.x = -2.0;
		topLeftFront = Point(x1, y1, z1);
		bottomRightBack = Point(x2, y2, z2);
		for (int i = 0; i < 8; ++i)
		{
			children[i] = nullptr;
		}
		for (int i = TopLeftFront;
			 i <= BottomLeftBack;
			 ++i)
			*children[i] = OctreeNodeCPU();
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
		if (children[pos]->point.x == -2.0)
		{
			return children[pos]->find(x, y, z);
		}
		// If an empty node is encountered
		if (children[pos]->point.x == -1.0)
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
class Octree_CPU : public OctreeBase
{
public:
	redwood::lang::cpu::CpuDevice *_device;
	OctreeNodeCPU *_root;
	std::vector<redwood::lang::DeviceAllocation *> buffers;

	Octree_CPU(float x1, float y1, float z1, float x2, float y2, float z2, redwood::lang::cpu::CpuDevice *device)
	{
		std::cout << "allocate" << std::endl;
		_device = device;
		redwood::lang::DeviceAllocation *root_buf = new redwood::lang::DeviceAllocation();
		_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCPU)), true, true}, root_buf);
		void *mappedData;
		buffers.push_back(root_buf);
		root_buf->device->map(*root_buf, &mappedData);
		_root = static_cast<OctreeNodeCPU *>(mappedData);
		_root->handle = root_buf;
	}

	~Octree_CPU()
	{
		for (int i = 0; i < buffers.size(); ++i)
		{
			_device->dealloc_memory(*buffers[i]);
		}
	}

	void insert(float x, float y, float z)
	{
		std::cout << "insert" << std::endl;
		// If the point already exists in the octree
		if (find(x, y, z) == true)
		{
			return;
		}
		insert(_root, x, y, z);
	}
	bool find(float x, float y, float z)
	{
		return _root->find(x, y, z);
	}
	void insert(OctreeNodeCPU *node, float x, float y, float z)
	{
		// If the point is out of bounds
		if (x < node->topLeftFront.x || x > node->bottomRightBack.x || y < node->topLeftFront.y || y > node->bottomRightBack.y || z < node->topLeftFront.z || z > node->bottomRightBack.z)
		{
			return;
		}

		// Binary search to insert the point
		float midx = (node->topLeftFront.x + node->bottomRightBack.x) / 2;
		float midy = (node->topLeftFront.y + node->bottomRightBack.y) / 2;
		float midz = (node->topLeftFront.z + node->bottomRightBack.z) / 2;

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
		if (node->children[pos]->point.x == -2.0)
		{
			// std::cout<<"internal point"<<std::endl;
			insert(node->children[pos], x, y, z);
			return;
		}

		// If an empty node is encountered
		else if (node->children[pos]->point.x == -1)
		{
			// std::cout <<"empty point"<<std::endl;
			//_device->dealloc_memory(*node->children[pos]->handle);
			redwood::lang::DeviceAllocation *child_buf = new redwood::lang::DeviceAllocation();
			_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCPU)), true, true}, child_buf);
			void *mappedData;
			buffers.push_back(child_buf);
			child_buf->device->map(*child_buf, &mappedData);
			OctreeNodeCPU *node = static_cast<OctreeNodeCPU *>(mappedData);
			*node = OctreeNodeCPU(x, y, z);
			node->handle = child_buf;
			node->children[pos] = node;
			return;
		}
		else
		{
			float x_ = node->children[pos]->point.x,
				  y_ = node->children[pos]->point.y,
				  z_ = node->children[pos]->point.z;
			//_device->dealloc_memory(*node->children[pos]->handle);
			node->children[pos] = nullptr;
			if (pos == TopLeftFront)
			{
				redwood::lang::DeviceAllocation *child_buf = new redwood::lang::DeviceAllocation();
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCPU)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				child_buf->device->map(*child_buf, &mappedData);
				OctreeNodeCPU *node = static_cast<OctreeNodeCPU *>(mappedData);
				*node = OctreeNodeCPU(node->topLeftFront.x,
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
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCPU)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				child_buf->device->map(*child_buf, &mappedData);
				OctreeNodeCPU *node = static_cast<OctreeNodeCPU *>(mappedData);
				*node = OctreeNodeCPU(midx + 1,
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
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCPU)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				child_buf->device->map(*child_buf, &mappedData);
				OctreeNodeCPU *node = static_cast<OctreeNodeCPU *>(mappedData);
				*node = OctreeNodeCPU(midx + 1,
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
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCPU)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				child_buf->device->map(*child_buf, &mappedData);
				OctreeNodeCPU *node = static_cast<OctreeNodeCPU *>(mappedData);
				*node = OctreeNodeCPU(node->topLeftFront.x,
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
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCPU)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				child_buf->device->map(*child_buf, &mappedData);
				OctreeNodeCPU *node = static_cast<OctreeNodeCPU *>(mappedData);
				*node = OctreeNodeCPU(node->topLeftFront.x,
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
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCPU)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				child_buf->device->map(*child_buf, &mappedData);
				OctreeNodeCPU *node = static_cast<OctreeNodeCPU *>(mappedData);
				*node = OctreeNodeCPU(midx + 1,
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
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCPU)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				child_buf->device->map(*child_buf, &mappedData);
				OctreeNodeCPU *node = static_cast<OctreeNodeCPU *>(mappedData);
				*node = OctreeNodeCPU(midx + 1,
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
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCPU)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				child_buf->device->map(*child_buf, &mappedData);
				OctreeNodeCPU *node = static_cast<OctreeNodeCPU *>(mappedData);
				*node = OctreeNodeCPU(node->topLeftFront.x,
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