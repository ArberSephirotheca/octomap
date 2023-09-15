#ifndef OCCUPANCY_MAP_H
#define OCCUPANCY_MAP_H

#include "octree.hpp"
#include "morton_util.hpp"
#include <librealsense2/rs.hpp>
#include <librealsense2/types.h>
using Code_t = uint64_t;

namespace map
{
    template <typename DATA_TYPE>
    class OccupancyMap : public oct::Octree<DATA_TYPE>
    {
    public:
        // constructor
        OccupancyMap(float resolution, int depth, float occupied_thres, float free_thres, float prob_hit, float prob_miss);
        // destructor
        virtual ~OccupancyMap() {}

        void insertPointCloud(const rs2::vertex &sensor_origin, rs2::points cloud, float max_range = -1,
                              int depth = 8)
        {
            std::vector<std::pair<Code_t, float>> occupied_hits;
            occupied_hits.reserve(cloud.size());
            rs2::points discretized;
            const rs2::vertex *points = cloud.get_vertices();
            int num_points = sizeof(points) / sizeof(points[0]);
            for (int i = 0; i < num_points; ++i)
            {
                rs::vertex end = points[i];
                rs2::vertex origin = sensor_origin;
                rs2::vertex direction = end - origin;
                float distance = direction.norm();

                if (moveLineInside(origin, end))
                {
                    continue;
                }

                if (0 > max_range || distance <= max_range)
                {
                    // Occupied space
                    Code_t end_code = PointToCode(end);
                    // TODO: indices to speedup
                    occupied_hits.push_back(std::make_pair(end_code, prob_hit_log_));
                }
                else
                {
                    direction /= distance;
                    end = origin + (direction * max_range);
                }

                discretized_push_back(end);

                for (int i : {0, 1, 2})
                {
                    min_change[i] = std::min(min_change[i], std::min(end[i], origin[i]));
                    max_change[i] = std::max(max_change[i], std::max(end[i], origin[i]));
                }
            }

            float prb_miss_log = prob_miss_log_ / float((2.0 * depth) + 1);
            // indices_.clear();

            insertPointCloudHelper(sensor_origin, std::move(discretized), 
                                    std::move(occupied_hits), prob_miss_log, 
                                    depth, min_change, max_change);
        }
        
        void insertPointCloudHelper(rs2::vertex snesor_origin, rs2::points&& discretized,
                                    std::vector<std::pair<Code_t, float>>&& occupied_hits,
                                    float prob_miss_log, int depth,
                                    rs2::vertex min_change, rs2::vertex max_change)
        {
            std::future<void> f = std::async(std::launch::async, [this, &occupied_hits](){
                std::for_each((occupied_hits.begin()), occupied_hits.end(),
                [this](auto&& hit){updateValue(hit.first, hit.second);});
            })
            CodeMap free_hits;
            freeSpace(sensor_origin, discretized, fre_hits, prob_miss_log, depth, free_hits);

            f.wait();

            for(const auto& [code, value] : free_hits){
                updateValue(code, value);
            }
        }

        void freeSpace(rs2::vertex snesor_origin, rs2::points&& discretized, float prob_miss_log, int depth, CodeMap& free_hits){

        }

        void updateValue(Code_t code, float ){

        }


    protected:
        float occupied_thres_;
        float free_thres_;
        float prob;
        // std::vector<Code_t> indices_;
    }
}

#endif // OCCUPANCY_MAP_H