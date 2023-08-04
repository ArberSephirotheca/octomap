#include <iostream>
#include <random>
#include "OctreeBase.cpp"

int main(){

    std::mt19937 generator(time(nullptr)); 
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    Octomap::Octree<Point> octree = Octomap::Octree<Point>();

    std::array<Point, 10000> pointcloud;
    for(auto& point : pointcloud){
        point.x = distribution(generator);
        point.y = distribution(generator);
        point.z = distribution(generator);
    }

    Code encoder = Code();
    std::array<OctreeLeafNode<Point>*, 16> path;
    for(auto point : pointcloud){
        octree.insertNode(encoder.morton3D(point), path, 16);
    }
}