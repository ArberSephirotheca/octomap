#ifndef OCCUPANCY_MAP_H
#define OCCUPANCY_MAP_H

#include "octree.hpp"
namespace map{
    class OccupancyMap : public oct::Octree{
        public:
            // constructor
            OccupancyMap(double resolution, int depth, double occupied_thres, double free_thres, double prob_hit, double prob_miss);
            // destructor
            virtual ~OccupancyMap(){}
        protected:``
            double occupied_thres_;
            double free_thres_;
            double prob;
    }
}

#endif // OCCUPANCY_MAP_H