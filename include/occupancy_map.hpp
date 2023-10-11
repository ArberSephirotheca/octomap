#ifndef OCCUPANCY_MAP_H
#define OCCUPANCY_MAP_H

#include "morton_util.hpp"
#include "octree.hpp"
#include <librealsense2/rs.hpp>
#include <pcl/filters/passthrough.h>
#include <pcl/point_types.h>

using Code_t = uint64_t;
namespace map {
template <typename DATA_TYPE>
class OccupancyMap : public oct::Octree<DATA_TYPE> {
protected:
  using Base = oct::Octree<DATA_TYPE>;
  using CodeMap = std::unordered_map<Code_t, DATA_TYPE>;

public:
  // constructor
  OccupancyMap(float resolution, int depth, float occupied_thres,
               float free_thres, float prob_hit, float prob_miss);
  // destructor
  virtual ~OccupancyMap() {}

  void insertPointCloud(
      const pcl::PointXYZ &sensor_origin,
      std::vector<pcl::PointXYZ, Eigen::aligned_allocator<pcl::PointXYZ>>
          &cloud,
      float max_range = -1, int depth = 8) {
    std::vector<std::pair<Code_t, float>> occupied_hits;
    occupied_hits.reserve(cloud.size());
    rs2::points discretized;
    for (auto &end : cloud) {
      pcl::PointXYZ origin = sensor_origin;
      pcl::PointXYZ direction = end - origin;
      float distance = direction.norm();

      if (!Base::moveLineInside(origin, end)) {
        continue;
      }

      if (0 > max_range || distance <= max_range) {
        // Occupied space
        Code_t end_code = PointToCode(end.x, end.y, end.z);
        // TODO: indices to speedup
        occupied_hits.push_back(std::make_pair(end_code, prob_hit_log_));
      } else {
        direction /= distance;
        end = origin + (direction * max_range);
      }

      discretized_push_back(end);

      for (int i : {0, 1, 2}) {
        min_change[i] = std::min(min_change[i], std::min(end[i], origin[i]));
        max_change[i] = std::max(max_change[i], std::max(end[i], origin[i]));
      }
    }

    float prb_miss_log = prob_miss_log_ / float((2.0 * depth) + 1);
    // indices_.clear();

    insertPointCloudHelper(sensor_origin, std::move(discretized),
                           std::move(occupied_hits), prob_miss_log, depth,
                           min_change, max_change);
  }

  void
  insertPointCloudHelper(pcl::PointXYZ snesor_origin, pcl_ptr &&discretized,
                         std::vector<std::pair<Code_t, float>> &&occupied_hits,
                         float prob_miss_log, int depth,
                         pcl::PointXYZ min_change, pcl::PointXYZ max_change) {
    std::future<void> f = std::async(std::launch::async, [this,
                                                          &occupied_hits]() {
      std::for_each((occupied_hits.begin()), occupied_hits.end(),
                    [this](auto &&hit) { updateValue(hit.first, hit.second); });
    }) CodeMap free_hits;
    freeSpace(sensor_origin, discretized, fre_hits, prob_miss_log, depth,
              free_hits);

    f.wait();

    for (const auto &[code, value] : free_hits) {
      updateValue(code, value);
    }
  }

  template <typename T>
  void
  freeSpace(pcl::PointXYZ snesor_origin,
            std::vector<pcl::PointXYZ, Eigen::aligned_allocator<pcl::PointXYZ>>
                &&discretized,
            float prob_miss_log, const T &value, int depth,
            CodeMap &indices) const {
    for (const auto &point : discretized) {
      pcl::PointXYZ current = sensor_origin;
      pcl::PointXYZ end;

      if (!Base::moveLineInside(current, end)) {
        continue;
      }

      freeSpaceNormal(current, end, indices, value);
    }
  }

  template <typename T>
  void freeSpaceNormal(const pcl::PointXYZ &from, const pcl::PointXYZ &to,
                       CodeMap &indices, const T &value,
                       int depth = 0, ) const {
    // Do it backwards
    pcl::PointXYZ current = to;
    pcl::PointXYZ end = from;

    pcl::PointXYZ direction = end - current;
    double distance = direction.norm();
    direction /= distance;
    Code_t current_key;
    Code_t end_key;
    std::array<int, 3> step;
    pcl::PointXYZ t_delta;
    pcl::PointXYZ t_max;
    Base::computeRayInit(current, end, direction, current_key,
                         end_key, step, t_delta, t_max, depth);
                    
		if (current_key == end_key) {
			indices.try_emplace(PointToCode(current_key.x, current_key.y, current_key.z ), value);
			return;
		}

		unsigned int already_update_in_row = 0;
		do {
			if (indices.try_emplace(PointToCode(current_key.x, current_key.y, current_key.z), value).second) {
				already_update_in_row = 0;
			} else {
				++already_update_in_row;
				if (0 < early_stopping && already_update_in_row >= early_stopping) {
					break;
				}
			}
			Base::computeRayTakeStep(current_key, step, t_delta, t_max);
		} while (current_key != end_key && t_max.min() <= distance);
  }

  void updateValue(Code_t code, float value) {
    // TODO
  }

protected:
  float occupied_thres_;
  float free_thres_;
  float prob_hit_log_;
  float prob_miss_log_;
  // std::vector<Code_t> indices_;
};
} // namespace map

#endif // OCCUPANCY_MAP_H