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

	// Represent the boundary of the cube
	redwood::lang::DeviceAllocationId children[8];
	Point topLeftFront, bottomRightBack;
	std::vector<redwood::lang::DeviceAllocation *> buffers;
	OctreeNodeCPU()
	{
		point = Point();
	}
	OctreeNodeCPU(int x, int y, int z)
	{
		point = Point(x, y, z);
	}
	OctreeNodeCPU(int x1, int y1, int z1, int x2, int y2, int z2, std::vector<redwood::lang::DeviceAllocation*>& buffers, redwood::lang::cpu::CpuDevice *device)
	{	
		// This use to construct Octree
		// with boundaries defined
		if (x2 < x1 || y2 < y1 || z2 < z1)
		{
			return;
		}
		point.x = -2;
		topLeftFront = Point(x1, y1, z1);
		bottomRightBack = Point(x2, y2, z2);
		redwood::lang::DeviceAllocation *child_buf = new redwood::lang::DeviceAllocation();
		device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCPU)*8), true, true}, child_buf);
		void *mappedData;
		// copy the pointer to buffer to the host
		buffers.push_back(child_buf);
		// maps the allocated buffer to host
		device->map(*child_buf, &mappedData);   
		for (int i = 0; i < 8; ++i)
		{
			auto current_node = OctreeNodeCPU();
			//current_node.handle = child_buf;
			// store the id of the buffer the child node is stored in
			children[i] = buffers.size()-1;
			// memcpy the data into mapped buffer
			memcpy(mappedData+i*(sizeof(OctreeNodeCPU)), &current_node, sizeof(OctreeNodeCPU));
			device->unmap(*child_buf);
		}
	}
	bool find(int x, int y, int z, int& count, redwood::lang::cpu::CpuDevice* device, std::vector<redwood::lang::DeviceAllocation *>& buffers)
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

		//std::cout<<"buffer size: "<<buffers.size()<<std::endl;
		//std::cout <<"id: " << children[pos]<<std::endl;
		void *mappedData;
		device->map(*buffers[children[pos]], &mappedData);
		OctreeNodeCPU *child_node = static_cast<OctreeNodeCPU *>(mappedData+pos*(sizeof(OctreeNodeCPU)));
		// If an internal node is encountered
		if (child_node->point.x == -2)
		{
			//std::cout<<"internal point"<<std::endl;
			count += 1;
			return child_node->find(x, y, z, count, device, buffers);
		}
		// If an empty node is encountered
		if (child_node->point.x == -1)
		{	
			//device->unmap(*buffers[children[pos]]);
			return false;
		}
		else
		{
			// If node is found with
			// the given value
			if (x == child_node->point.x && y == child_node->point.y && z == child_node->point.z){
			// This unmap causes incorrect result
				//device->unmap(*buffers[children[pos]]);
				return 1;
			}
		}
		device->unmap(*buffers[children[pos]]);
		return 0;
	}
};

// Octree class
class Octree_CPU : public OctreeBase
{
public:
	redwood::lang::cpu::CpuDevice *_device;
	redwood::lang::DeviceAllocation *_root;
	std::vector<redwood::lang::DeviceAllocation *> buffers;
	int find_count = 0;

	Octree_CPU(int x1, int y1, int z1, int x2, int y2, int z2, redwood::lang::cpu::CpuDevice *device)
	{
		// TODO: This is a bug
		//std::cout << "allocate" << std::endl;
		_device = device;
		redwood::lang::DeviceAllocation *root_buf = new redwood::lang::DeviceAllocation();
		_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCPU)), true, true}, root_buf);
		void *mappedData;
		buffers.push_back(root_buf);
		root_buf->device->map(*root_buf, &mappedData);
		auto root_node = OctreeNodeCPU(x1, y1, z1, x2, y2, z2, buffers, _device);
		//root_node.handle = root_buf;
		memcpy(mappedData, &root_node, sizeof(OctreeNodeCPU));
		//_device->unmap(*root_buf);
		_root = root_buf;
		//std::cout<<"done allocate"<<std::endl;
	}

	~Octree_CPU()
	{
		while (!buffers.empty())
		{
			//std::cout<<"buffer size: "<<buffers.size()<<std::endl;
			_device->dealloc_memory(*buffers.back());
			buffers.pop_back();
		}
	}

	void insert(int x, int y, int z)
	{
		// If the point already exists in the octree
		if (find(x, y, z) == true)
		{
			//std::cout<<"already exist"<<std::endl;
			return;
		}
		void *currentData;
		_device->map(*_root, &currentData);
		OctreeNodeCPU *root_node = static_cast<OctreeNodeCPU *>(currentData);
		insert(root_node, x, y, z);
		//_device->unmap(*_root);
	}
	bool find(int x, int y, int z)
	{
		//std::cout<<"find"<<std::endl;
		int count = 0;
		void *currentData;
		_device->map(*_root, &currentData);
		OctreeNodeCPU *root_node = static_cast<OctreeNodeCPU *>(currentData);
		bool result =  root_node->find(x, y, z,count, _device, buffers);
		find_count += count;
		//_device->unmap(*_root);
		return result;
	}
	void insert(OctreeNodeCPU* root_node, int x, int y, int z)
	{
		//std::cout<<"insert"<<std::endl;
		// If the point is out of bounds
		if (x < root_node->topLeftFront.x || x > root_node->bottomRightBack.x || y < root_node->topLeftFront.y || y > root_node->bottomRightBack.y || z < root_node->topLeftFront.z || z > root_node->bottomRightBack.z)
		{
			return;
		}

		// Binary search to insert the point
		int midx = (root_node->topLeftFront.x + root_node->bottomRightBack.x) / 2;
		int midy = (root_node->topLeftFront.y + root_node->bottomRightBack.y) / 2;
		int midz = (root_node->topLeftFront.z + root_node->bottomRightBack.z) / 2;

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
		// maps child node to the host
		void *childData;
		_device->map(*buffers[root_node->children[pos]], &childData);
		OctreeNodeCPU *child_node = static_cast<OctreeNodeCPU *>(childData+pos*(sizeof(OctreeNodeCPU)));

		// If an internal node is encountered
		if (child_node->point.x == -2)
		{
			// std::cout<<"internal point"<<std::endl;
			insert(child_node, x, y, z);
			//_device->unmap(*buffers[root_node->children[pos]]);
			return;
		}

		// If an empty node is encountered
		else if (child_node->point.x == -1)
		{
			 //std::cout <<"empty point"<<std::endl;
			 // allocate memory for children node
			//_device->dealloc_memory(*buffers[root_node->children[pos]]);
			redwood::lang::DeviceAllocation *child_buf = new redwood::lang::DeviceAllocation();
			_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCPU)), true, true}, child_buf);
			void *mappedData;
			// copy the pointer to buffer to the host
			buffers.push_back(child_buf);
			// maps the allocated buffer to host
			root_node->children[pos] = buffers.size()-1;
			_device->map(*child_buf, &mappedData);
			auto children_node = OctreeNodeCPU(x, y, z);
			//children_node.handle = child_buf;
			// store the id of the buffer the child node is stored in
			// memcpy the data into mapped buffer
			memcpy(mappedData, &children_node, sizeof(OctreeNodeCPU));
			// unmap the buffer
			//_device->unmap(*child_buf);
			//_device->unmap(*buffers[root_node->children[pos]]);
			return;
		}
		else
		{
			//std::cout<<"point"<<std::endl;
			int x_ = child_node->point.x,
				  y_ = child_node->point.y,
				  z_ = child_node->point.z;
			//_device->dealloc_memory(*child_node->handle);
			child_node = nullptr;
			//_device->unmap(*buffers[root_node->children[pos]]);
			if (pos == TopLeftFront)
			{
				redwood::lang::DeviceAllocation *child_buf = new redwood::lang::DeviceAllocation();
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCPU)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				root_node->children[pos] = buffers.size()-1;
				_device->map(*child_buf, &mappedData);
				auto children_node = OctreeNodeCPU(root_node->topLeftFront.x,
									  root_node->topLeftFront.y,
									  root_node->topLeftFront.z,
									  midx,
									  midy,
									  midz,  buffers, _device);
				//children_node.handle = child_buf;
				memcpy(mappedData, &children_node, sizeof(OctreeNodeCPU));
				child_node = static_cast<OctreeNodeCPU *>(mappedData);
				//_device->unmap(*child_buf);
			}

			else if (pos == TopRightFront)
			{
				redwood::lang::DeviceAllocation *child_buf = new redwood::lang::DeviceAllocation();
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCPU)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				root_node->children[pos] = buffers.size()-1;
				_device->map(*child_buf, &mappedData);
				auto children_node = OctreeNodeCPU(midx + 1,
									  root_node->topLeftFront.y,
									  root_node->topLeftFront.z,
									  root_node->bottomRightBack.x,
									  midy,
									  midz,  buffers, _device);
				//children_node.handle = child_buf;
				memcpy(mappedData,  &children_node, sizeof(OctreeNodeCPU));
				child_node = static_cast<OctreeNodeCPU *>(mappedData);
				//_device->unmap(*child_buf);

			}
			else if (pos == BottomRightFront)
			{
				redwood::lang::DeviceAllocation *child_buf = new redwood::lang::DeviceAllocation();
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCPU)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				root_node->children[pos] = buffers.size()-1;
				_device->map(*child_buf, &mappedData);
				auto children_node = OctreeNodeCPU(midx + 1,
									  midy + 1,
									  root_node->topLeftFront.z,
									  root_node->bottomRightBack.x,
									  root_node->bottomRightBack.y,
									  midz,  buffers, _device);
				//children_node.handle = child_buf;
				memcpy(mappedData,  &children_node, sizeof(OctreeNodeCPU));
				child_node = static_cast<OctreeNodeCPU *>(mappedData);
				//_device->unmap(*child_buf);
			}
			else if (pos == BottomLeftFront)
			{
				redwood::lang::DeviceAllocation *child_buf = new redwood::lang::DeviceAllocation();
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCPU)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				root_node->children[pos] = buffers.size()-1;
				_device->map(*child_buf, &mappedData);
				auto children_node = OctreeNodeCPU(root_node->topLeftFront.x,
									  midy + 1,
									  root_node->topLeftFront.z,
									  midx,
									  root_node->bottomRightBack.y,
									  midz,  buffers, _device);
				//children_node.handle = child_buf;
				memcpy(mappedData,  &children_node, sizeof(OctreeNodeCPU));
				child_node = static_cast<OctreeNodeCPU *>(mappedData);
				//_device->unmap(*child_buf);
			}
			else if (pos == TopLeftBottom)
			{
				redwood::lang::DeviceAllocation *child_buf = new redwood::lang::DeviceAllocation();
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCPU)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				root_node->children[pos] = buffers.size()-1;
				_device->map(*child_buf, &mappedData);
				auto children_node = OctreeNodeCPU(root_node->topLeftFront.x,
									  root_node->topLeftFront.y,
									  midz + 1,
									  midx,
									  midy,
									  root_node->bottomRightBack.z,  buffers, _device);
				//children_node.handle = child_buf;
				memcpy(mappedData,  &children_node, sizeof(OctreeNodeCPU));
				child_node = static_cast<OctreeNodeCPU *>(mappedData);
				//_device->unmap(*child_buf);
			}
			else if (pos == TopRightBottom)
			{
				redwood::lang::DeviceAllocation *child_buf = new redwood::lang::DeviceAllocation();
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCPU)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				root_node->children[pos] = buffers.size()-1;
				_device->map(*child_buf, &mappedData);
				auto children_node = OctreeNodeCPU(midx + 1,
									  root_node->topLeftFront.y,
									  midz + 1,
									  root_node->bottomRightBack.x,
									  midy,
									  root_node->bottomRightBack.z,  buffers, _device);
				//children_node.handle = child_buf;
				memcpy(mappedData,  &children_node, sizeof(OctreeNodeCPU));
				child_node = static_cast<OctreeNodeCPU *>(mappedData);

				//_device->unmap(*child_buf);
			}
			else if (pos == BottomRightBack)
			{
				redwood::lang::DeviceAllocation *child_buf = new redwood::lang::DeviceAllocation();
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCPU)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				root_node->children[pos] = buffers.size()-1;
				_device->map(*child_buf, &mappedData);
				auto children_node = OctreeNodeCPU(midx + 1,
									  midy + 1,
									  midz + 1,
									  root_node->bottomRightBack.x,
									  root_node->bottomRightBack.y,
									  root_node->bottomRightBack.z,  buffers, _device);
				//children_node.handle = child_buf;
				memcpy(mappedData,  &children_node, sizeof(OctreeNodeCPU));
				child_node = static_cast<OctreeNodeCPU *>(mappedData);
				//_device->unmap(*child_buf);

			}
			else if (pos == BottomLeftBack)
			{
				redwood::lang::DeviceAllocation *child_buf = new redwood::lang::DeviceAllocation();
				_device->allocate_memory(redwood::lang::Device::AllocParams{(sizeof(OctreeNodeCPU)), true, true}, child_buf);
				void *mappedData;
				buffers.push_back(child_buf);
				root_node->children[pos] = buffers.size()-1;
				_device->map(*child_buf, &mappedData);
				auto children_node = OctreeNodeCPU(root_node->topLeftFront.x,
									  midy + 1,
									  midz + 1,
									  midx,
									  root_node->bottomRightBack.y,
									  root_node->bottomRightBack.z,  buffers, _device);
				//children_node.handle = child_buf;
				memcpy(mappedData,  &children_node, sizeof(OctreeNodeCPU));
				child_node = static_cast<OctreeNodeCPU *>(mappedData);
				//_device->unmap(*child_buf);
			}
			insert(child_node, x_, y_, z_);
			insert(child_node, x, y, z);
			 //_device->unmap(*buffers[root_node->children[pos]]);
		}
	}
};