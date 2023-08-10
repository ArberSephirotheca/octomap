#include "OctreeBase.cpp"
#include <chrono>
#include <cstdlib>
#include <iostream>
int main() {

  Octomap::Octree<Point> octree = Octomap::Octree<Point>();

  std::vector<Point> pointcloud;
  for (int i = 0; i < 280000; i++) {
    Point point = Point{};
    point.x = static_cast<double>(rand()) / RAND_MAX;
    point.y = static_cast<double>(rand()) / RAND_MAX;
    point.z = static_cast<double>(rand()) / RAND_MAX;
    pointcloud.push_back(point);
  }

  Code encoder = Code();
  auto start = std::chrono::high_resolution_clock::now();
  for (auto point : pointcloud) {
    Code *code = new Code(encoder.morton3D(point), 0);

    octree.insertNode(*code);
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count();

  std::cout << "Time taken for insertion: " << duration << " milliseconds"
            << std::endl;

  start = std::chrono::high_resolution_clock::now();
  for (auto point : pointcloud) {
    Code code = Code(encoder.morton3D(point), 0);
    // std::cout <<"search code: "<< code.code_<<std::endl;
    // std::cout <<"search depth: "<< code.depth_<<std::endl;
    const Point *result = octree.search(code);
  }

  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                 .count();

  std::cout << "Time taken for search: " << duration << " milliseconds"
            << std::endl;

  std::cout << "numer of leaf nodes: " << octree.getNumLeafNodes() << std::endl;
  std::cout << "numer of inner nodes: " << octree.getNumInnerNodes()
            << std::endl;
  std::cout << "numer of inner leaf nodes: " << octree.getNumInnerLeafNodes()
            << std::endl;
}