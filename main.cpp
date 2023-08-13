#include "OctreeBase.cpp"
#include "binary_radix_tree.cpp"
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <bitset>

    void test_basic(){
    
        std::string binaryString = "0000000000000000000000000011100";
    std::bitset<32> binaryNumber(binaryString);
    
    uint32_t result = binaryNumber.to_ulong();
    std::cout << "result: " <<result<<std::endl;
    std::cout<<"30: "<<std::bitset<32>(30)<<std::endl;
		const std::vector<uint32_t> morton_keys{ 1, 2, 4, 5, 19, 24, 25, 28, 29, 30}; // the 8 keys of the same value as in the paper
		const auto brt{ binary_radix_tree::create(10, morton_keys, 8)};

		const bool allright{	brt.internal_nodes[0] == 3 && brt.internal_nodes[1] == 4 &&
								brt.internal_nodes[6] == 1 && brt.internal_nodes[7] == 2 &&
								brt.internal_nodes[9] == 5 && brt.internal_nodes[10] == 6 };		// expected values
    for(int i = 0; i < 30; i++){
      std::cout <<"nodes[" <<i<<"]"<<brt.internal_nodes[i]<<std::endl;
      std::cout <<"bitset: "<<   std::bitset<32>(brt.internal_nodes[i])<<std::endl;
    }
		std::cout << "paper test " << ((allright) ? "ok." : "not ot.") << std::endl;
    };
int main() {

    test_basic();

  Octomap::Octree<Point> octree = Octomap::Octree<Point>();

  std::vector<Point> pointcloud{};
  std::vector<uint32_t> sortedMorton{};
  for (int i = 0; i < 1; ++i) {
    Point point = Point{};
    point.x = static_cast<double>(rand()) / RAND_MAX;
    point.y = static_cast<double>(rand()) / RAND_MAX;
    point.z = static_cast<double>(rand()) / RAND_MAX;
    pointcloud.push_back(point);
  }

  Code encoder = Code();
  std::cout <<"size of pointcloud: "<< pointcloud.size()<<std::endl;
  auto start = std::chrono::high_resolution_clock::now();
  for (const Point& point : pointcloud) {
    uint32_t mortonCode = encoder.morton3D(point);
    std::cout <<"point x:" << point.x<<std::endl;
        std::cout <<"point y:" << point.x<<std::endl;
            std::cout <<"point z:" << point.x<<std::endl;
    std::cout <<"morton code: "<< std::bitset<32>(mortonCode)<<std::endl;
    Code* code = new Code(mortonCode, 0);
    octree.insertNode(*code);
    sortedMorton.push_back(mortonCode);
  }
  radixSortMorton(sortedMorton);

  binary_radix_tree::create(sortedMorton.size(),sortedMorton, 8);

  

  
  auto end = std::chrono::high_resolution_clock::now();
  auto duration =std::chrono::duration_cast<std::chrono::milliseconds>(end - start)    .count();

  std::cout << "Time taken for insertion: " << duration << " milliseconds"     << std::endl;

  start = std::chrono::high_resolution_clock::now();
  for (auto point : pointcloud) {
    Code code = Code(encoder.morton3D(point), 0);
     std::cout <<"search code: "<< std::bitset<32>(code.getCode())<<std::endl;
     std::cout <<"search depth: "<< code.getDepth()<<std::endl;
    const Point *result = octree.search(code);
  }

  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start)              .count();

  std::cout << "Time taken for search: " << duration << " milliseconds"         << std::endl;

  std::cout << "numer of leaf nodes: " << octree.getNumLeafNodes() << std::endl;
  std::cout << "numer of inner nodes: " << octree.getNumInnerNodes()
            << std::endl;
  std::cout << "numer of inner leaf nodes: " << octree.getNumInnerLeafNodes()
            << std::endl;
            
}