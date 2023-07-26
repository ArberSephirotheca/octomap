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
        Octree(int map_scale[2], double voxel_scale, int min_occup_threshold, double min_ray_length, double max_ray_length, int max_submap_num, int K, double min_x, double min_y,
               double min_z, double max_x, double max_y, double max_z)
        {
            this->max_size_xy = map_scale[0];
            this->max_size_z = map_scale[1];
            this->voxel_scale = voxel_scale;
            this->min_occup_threshold = min_occup_threshold;
            this->min_ray_length = min_ray_length;
            this->max_ray_length = max_ray_length;
            this->max_submap_num = max_submap_num;
            this->K = K;
            std::unique_ptr<ChunkedAllocator<uint32_t>> octreeAllocator(new ChunkedAllocator<uint32_t>());
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
            if (!containsPoint(minBound, maxBound, point))
            {
                return; // Point is outside the octree bounds
            }
            insertPoint(root, minBound, maxBound, point, 0);
        }

        void query(const Point &queryPoint, std::vector<Point> &result)
        {
            if (!containsPoint(minBound, maxBound, queryPoint))
            {
                return; // Query point is outside the octree bounds
            }
            queryPointInNode(root, minBound, maxBound, queryPoint, result);
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
        std::vector<OctreeNode> nodes;
        OctreeNode *root;

        bool containsPoint(const Point &minPoint, const Point &maxPoint, const Point &point) const
        {
            return point.x >= minPoint.x && point.y >= minPoint.y && point.z >= minPoint.z &&
                   point.x <= maxPoint.x && point.y <= maxPoint.y && point.z <= maxPoint.z;
        }

        OctreeNode *createNode(const Point &minPoint, const Point &maxPoint, int depth)
        {
            nodes.emplace_back();
            if (depth == MAX_DEPTH)
            {
                return &nodes.back();
            }
            
            nodes.back().points.reserve(MAX_OBJECTS_PER_NODE);
            double midX = (minPoint.x + maxPoint.x) / 2.0;
            double midY = (minPoint.y + maxPoint.y) / 2.0;
            double midZ = (minPoint.z + maxPoint.z) / 2.0;

            int octantIndex = ((mortonEncode(minPoint) >= mortonEncode(Point(midX, midY, midZ))) << 2) |
                              ((mortonEncode(minPoint) >= mortonEncode(Point(midX, midY, maxPoint.z))) << 1) |
                              (mortonEncode(minPoint) >= mortonEncode(Point(midX, maxPoint.y, midZ)));
            for (int i = 0; i < 8; ++i)
            {
                if (i == octantIndex)
                {
                    nodes.back().children[i] = nodes.size() - 1; // Set the child index for the current octant
                }
                else
                {
                    nodes.back().children[i] = -1; // Initialize other children to -1 (invalid index)
                }
            }
            return &nodes.back();
        }
        // Function to compute the Morton code (Z-order curve) for a 3D point (x, y, z)
        unsigned int mortonEncode(const Point &point)
        {
            unsigned int result = 0;
            unsigned x = point.x;
            unsigned y = point.y;
            unsigned z = point.z;
            for (unsigned int i = 0; i < 10; ++i)
            {
                result |= ((x & (1u << i)) << 2 * i) | ((y & (1u << i)) << (2 * i + 1)) | ((z & (1u << i)) << (2 * i + 2));
            }
            return result;
        }

        // Function to decode the 3D coordinates (x, y, z) from a Morton code
        void mortonDecode(unsigned int morton, const Point &point)
        {
            unsigned x = point.x;
            unsigned y = point.y;
            unsigned z = point.z;
            for (unsigned int i = 0; i < 10; ++i)
            {
                x |= ((morton & (1u << (3 * i))) >> (2 * i));
                y |= ((morton & (1u << (3 * i + 1))) >> (2 * i + 1));
                z |= ((morton & (1u << (3 * i + 2))) >> (2 * i + 2));
            }
        }
        Point octantToMinPoint(const Point &minPoint, const Point &maxPoint, int octantIndex) const
        {
            double midX = (minPoint.x + maxPoint.x) / 2.0;
            double midY = (minPoint.y + maxPoint.y) / 2.0;
            double midZ = (minPoint.z + maxPoint.z) / 2.0;

            double minX = (octantIndex & 1) ? midX : minPoint.x;
            double minY = (octantIndex & 2) ? midY : minPoint.y;
            double minZ = (octantIndex & 4) ? midZ : minPoint.z;
            return Point(minX, minY, minZ);
        }

        Point octantToMaxPoint(const Point &minPoint, const Point &maxPoint, int octantIndex) const
        {
            double midX = (minPoint.x + maxPoint.x) / 2.0;
            double midY = (minPoint.y + maxPoint.y) / 2.0;
            double midZ = (minPoint.z + maxPoint.z) / 2.0;

            double maxX = (octantIndex & 1) ? maxPoint.x : midX;
            double maxY = (octantIndex & 2) ? maxPoint.y : midY;
            double maxZ = (octantIndex & 4) ? maxPoint.z : midZ;
            return Point(maxX, maxY, maxZ);
        }
        void insertPoint(OctreeNode *node, const Point &minPoint, const Point &maxPoint, const Point &point, int depth)
        {
            if (node->points.size() < MAX_OBJECTS_PER_NODE || depth == MAX_DEPTH)
            {
                node->points.push_back(point);
                return;
            }

            double midX = (minPoint.x + maxPoint.x) / 2.0;
            double midY = (minPoint.y + maxPoint.y) / 2.0;
            double midZ = (minPoint.z + maxPoint.z) / 2.0;

            int octantIndex = ((mortonEncode(point) >= mortonEncode(Point(midX, midY, midZ))) << 2) |
                              ((mortonEncode(point) >= mortonEncode(Point(midX, midY, maxPoint.z))) << 1) |
                              (mortonEncode(point) >= mortonEncode(Point(midX, maxPoint.y, midZ)));

            if (node->children[octantIndex] == -1)
            {
                node->children[octantIndex] = nodes.size(); // Set the child index for the current octant
                nodes.push_back(OctreeNode());              // Add a new node for the current octant
                nodes.back().points.reserve(MAX_OBJECTS_PER_NODE);
            }

            insertPoint(&nodes[node->children[octantIndex]], octantToMinPoint(minPoint, maxPoint, octantIndex),
                        octantToMaxPoint(minPoint, maxPoint, octantIndex), point, depth + 1);
        }

        void queryPointInNode(OctreeNode *node, const Point &minPoint, const Point &maxPoint, const Point &queryPoint, std::vector<Point> &result)
        {
            if (node == nullptr)
            {
                return;
            }

            if (node->points.size() > 0)
            {
                // Add points in the node to the result vector
                result.insert(result.end(), node->points.begin(), node->points.end());
            }

            double midX = (minPoint.x + maxPoint.x) / 2.0;
            double midY = (minPoint.y + maxPoint.y) / 2.0;
            double midZ = (minPoint.z + maxPoint.z) / 2.0;

            int octantIndex = ((mortonEncode(queryPoint) >= mortonEncode(Point(midX, midY, midZ))) << 2) |
                              ((mortonEncode(queryPoint) >= mortonEncode(Point(midX, midY, maxPoint.z))) << 1) |
                              (mortonEncode(queryPoint) >= mortonEncode(Point(midX, maxPoint.y, midZ)));

            queryPointInNode(&nodes[node->children[octantIndex]], octantToMinPoint(minPoint, maxPoint, octantIndex),
                             octantToMaxPoint(minPoint, maxPoint, octantIndex), queryPoint, result);
        }
    }
}