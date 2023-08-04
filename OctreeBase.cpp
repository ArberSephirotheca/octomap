#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <vector>
#include <stdint.h>
#include "ChunkedAllocator.h"
#include "Point.h"
#include "OctreeNode.h"
#include "code.h"

namespace Octomap
{
    const int MAX_POINTS_PER_NODE = 4;
    const int MAX_DEPTH = 16;
    const int MAX_OBJECTS_PER_NODE = 4;
    struct Block
    {
        std::array<std::array<double, 2>, 2> pixels;
    };

    template <typename DATA_TYPE, typename INNER_NODE = OctreeInnerNode<DATA_TYPE>,
              typename LEAF_NODE = OctreeLeafNode<DATA_TYPE>>
    class Octree
    {
    protected:
        using Path = std::array<LEAF_NODE*, MAX_DEPTH>;

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
            // nodes.push_back(OctreeNode(0, false));
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

        // return corresponding node and depth given the morton code
        std::pair<LEAF_NODE const *, int> getNode(const Code &code) const
        {
            LEAF_NODE const *node = &getRoot();
            // depth 0 = leaf node
            for (int depth = depth_ - 1; depth > code.getDepth(); --depth)
            {
                INNER_NODE const &inner_node = static_cast<INNER_NODE const &>(*node);
                // if reach leaf node, return
                if (inner_node.is_leaf)
                {
                    return std::make_pair(node, depth + 1);
                }
                // otherwise, move to next level
                node = &getChild(inner_node, depth, code.getChildIndex(depth));
            }
            return std::make_pair(node, code.getDepth());
        }

        std::pair<Path, int> getNodePath(const Code& code){
            Path path;
            path[depth_] = static_cast<LEAF_NODE*>(&getRoot());
            int depth = depth_;
            for(; depth > code.getDepth(); -- depth){
                INNER_NODE& node = static_cast<INNER_NODE&>(*path[depth]);
                // if node is leaf node, no need to traverse down
                if(node.is_leaf){
                    break;
                }

                int child_depth = depth-1;
                path[child_depth] = static_cast<LEAF_NODE*>(&getChild(node, child_depth, code.getChildIndex(child_depth)));
            }

            return std::make_pair(path, depth);
        }

        void insertNode(const Code& code, Path& path, int depth)
        {   
            // depth 0 = leaf node
            for(; depth > code.getDepth(); --depth){
                INNER_NODE& node = static_cast<INNER_NODE&>(*(path[depth]));
                // TODO: change is_leaf to hasChildren
                if(node.is_leaf){
                    createChildren(node, depth);
                }
                int child_depth = depth -1;
                path[child_depth] = static_cast<LEAF_NODE*>(&getChild(node, child_depth,code.getChildIndex(child_depth)));
            }

        }

        DATA_TYPE const* search(Code const &code) const
        {
            auto [node, depth] = getNode(code);
            if(depth == code.getDepth()){
                return &(node->value);
            }else{
                return nullptr;
            }
        }

        LEAF_NODE &getChild(INNER_NODE const &inner_node, int child_depth, int child_idx)
        {
            if (child_depth == 0)
            {
                return getLeafChild(inner_node, child_idx);
            }
            else
            {
                return getInnerChild(inner_node, child_idx);
            }
        }

        LEAF_NODE &getLeafChild(INNER_NODE const &inner_node, int child_idx)
        {
            return getLeafChildren(inner_node)[child_idx];
        }

        INNER_NODE &getInnerChild(INNER_NODE const &inner_node, int child_idx)
        {
            return getInnerChildren(inner_node)[child_idx];
        }

        std::array<LEAF_NODE &, 8> getLeafChildren(INNER_NODE const &inner_node)
        {
            return static_cast<LEAF_NODE>(inner_node.children);
        }

        std::array<INNER_NODE&, 8>getInnerChildren(INNER_NODE const &inner_node, int child_idx)
        {
            return static_cast<INNER_NODE>(inner_node.children);
        }
        bool createChildren(INNER_NODE &node, int depth){
            // Cannot rcreate children for non-leaf node
            if (!node.is_leaf){
                return false;
            }
            
            // if pointer to children nodes is null, create children
            if(!node.children){
                // second lowest level, create leaf node children
                if(1 == depth){
                    node.children = new std::array<LEAF_NODE, 8>();
                    num_leaf_nodes += 8;
                    num_inner_leaf_node -= 1;
                }else{
                    node.children = new std::array<INNER_NODE, 8>();
                    num_inner_leaf_node += 7;
                }
                num_inner_node += 1;
            }
            // TODO

            node.is_leaf = false;
            return true;
        }

        INNER_NODE const& getRoot()const{
            return root_;
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
        // 10 cm voxel size.
        double resolution = 0.1;
        // Number of levels for the octree. Has to be [2, 21].
        // This determines the maximum volume that can be represented by the map
        // resolution * 2^(depth_levels) meter in each dimension.
        int depth_ = 16;
        Point minBound;
        Point maxBound;
        ChunkedAllocator<Point> *allocator;
        std::vector<std::pair<int, int>> data_points; //[morton code, value]
        //Path path; // represent arrays of nodes for each level
        INNER_NODE root_;
        int num_leaf_nodes = 0;
        int num_inner_leaf_node = 0;
        int num_inner_node = 0;

        /*
        std::vector<int> queryPoint(int node_idx, uint32_t code, int depth)
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
                    return queryPoint(node.children[child_idx], code, depth + 1);
                    // if not set, return empty result
                }
                else
                {
                    return {};
                }
            }
        }
        */
       /*
        void sorted_morton_code()
        {
            std::sort(nodes.begin(), nodes.end(), [](const OctreeNode &a, const OctreeNode &b)
                      { return a.code < b.code; });
        }
        */
    };
}