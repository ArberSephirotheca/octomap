#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <stdint.h>
#include "ChunkedAllocator.hpp"
#include "Point.hpp"
#include "OctreeNode.hpp"

namespace Octomap
{
    const int MAX_POINTS_PER_NODE = 4;
    const int MAX_DEPTH = 16;
    const int MAX_OBJECTS_PER_NODE = 4;
    struct Block
    {
        std::array<std::array<double, 2>, 2> pixels;
    };

    class Octree
    {
    public:
        Octree(/*int map_scale[2], double voxel_scale, int min_occup_threshold, double min_ray_length, double max_ray_length, int max_submap_num, int K, double min_x, double min_y,
               double min_z, double max_x, double max_y, double max_z*/
        )
        { /*
          this->max_size_xy = map_scale[0];
          this->max_size_z = map_scale[1];
          this->voxel_scale = voxel_scale;
          this->min_occup_threshold = min_occup_threshold;
          this->min_ray_length = min_ray_length;
          this->max_ray_length = max_ray_length;
          this->max_submap_num = max_submap_num;
          this->K = K;
          */
            nodes.push_back(OctreeNode(0, 0));
            this->allocator = new ChunkedAllocator<Point>();
            this->initialize_field();
            this->construct_octree();
            this->initialize_submap_fields();
        }

        void initialize_field()
        {
        }

        void construct_octree()
        {
            int childCount = 0;
            int childIndices[8];
        }

        void initialize_submap_fields()
        {
        }

        bool is_occupy()
        {
        }

        void recast_pcl_to_map()
        {
        }

        void recast_depth_to_map()
        {
        }

        void fuse_submaps()
        {
        }

        std::vector<double> get_occup_voxels()
        {
        }

        void insert(const Point &point)
        {
            uint32_t code = morton3D(point);
            insertPoint(0, code, point, 0);
        }

        std::vector<int> query(const Point &point)
        {
            uint32_t code = morton3D(point);
            return queryPoint(0, code, point, 0);
        }

    private:
        double **map;
        int max_size_xy;
        int max_size_z;
        double voxel_scale;
        int min_occup_threshold;
        double min_ray_length;
        double max_ray_length;
        int max_submap_num;
        int K;
        Point minBound;
        Point maxBound;
        ChunkedAllocator<Point> *allocator;
        std::vector<std::pair<int, int>> data_points;
        std::vector<OctreeNode> nodes;
        OctreeNode *root;
        
        // see https://developer.nvidia.com/blog/thinking-parallel-part-iii-tree-construction-gpu/
        unsigned int expandBits(unsigned int v)
        {
            v = (v * 0x00010001u) & 0xFF0000FFu;
            v = (v * 0x00000101u) & 0x0F00F00Fu;
            v = (v * 0x00000011u) & 0xC30C30C3u;
            v = (v * 0x00000005u) & 0x49249249u;
            return v;
        }

        // Calculates a 30-bit Morton code for the
        // given 3D point located within the unit cube [0,1].
        unsigned int morton3D(const Point& point)
        {   
            double x = point.x;
            double y = point.y;
            double z = point.z;
            x = std::min(std::max(x * 1024.0, 0.0), 1023.0);
            y = std::min(std::max(y * 1024.0, 0.0), 1023.0);
            z = std::min(std::max(z * 1024.0, 0.0), 1023.0);
            unsigned int xx = expandBits((unsigned int)x);
            unsigned int yy = expandBits((unsigned int)y);
            unsigned int zz = expandBits((unsigned int)z);
            return xx * 4 + yy * 2 + zz;
        }

        void insertPoint(int node_idx, uint32_t code, const Point &point, int depth)
        {
            OctreeNode &node = nodes[node_idx];
            if (depth == MAX_DEPTH)
            {
                if (node.data_index == -1)
                {
                    node.data_index = data_points.size();
                }
                node.data_index++;
                data_points.push_back(std::make_pair(1, depth));
            }
            else
            {
                // get child index by depth and morton code
                int child_idx = get_child_index(code, depth);
                std::cout<<"child idx:" << child_idx<<std::endl;
                // if child not set, initlaize the child and place last index on the entry.
                if (node.children[child_idx] == -1)
                {
                    node.children[child_idx] = nodes.size();
                    // calculate the child code
                    uint32_t child_code = code | (child_idx << (3 * (MAX_DEPTH - depth - 1)));
                    nodes.push_back(OctreeNode(child_code, depth + 1));
                }
                insertPoint(node.children[child_idx], code, point, depth + 1);
            }
        }

        std::vector<int> queryPoint(int node_idx, uint32_t code, const Point &point, int depth)
        {
            const OctreeNode &node = nodes[node_idx];
            if (depth == MAX_DEPTH)
            {
                std::vector<int> result;
                for (int i = 0; i < node.data_index; ++i)
                {
                    result.push_back(data_points[node.data_index + i].first);
                }
                return result;
            }
            else
            {
                int child_idx = get_child_index(code, depth);
                if (node.children[child_idx] != -1)
                {
                    return queryPoint(node.children[child_idx], code, point, depth + 1);
                    // if not set, return empty result
                }
                else
                {
                    return {};
                }
            }
        }

        int get_child_index(uint32_t code, int depth)
        {
            return (code >> (3 * (MAX_DEPTH - depth - 1))) & 0x7;
        }
    };
}