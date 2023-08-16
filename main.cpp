#include "OctreeBase.cpp"
#include "binary_radix_tree.cpp"
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <bitset>

    void test_basic(std::vector<uint32_t>& morton_keys){
    radixSortMorton(morton_keys);
    int key_num = morton_keys.size();
		//const std::vector<uint32_t> morton_keys{ 1, 2, 4, 5, 19, 24, 25, 30}; // the 8 keys of the same value as in the paper
		const auto brt{ binary_radix_tree::create(key_num, morton_keys, 8)};

		const bool allright{	brt.internal_nodes[0] == 3 && brt.internal_nodes[1] == 4 &&
								brt.internal_nodes[6] == 1 && brt.internal_nodes[7] == 2 &&
								brt.internal_nodes[9] == 5 && brt.internal_nodes[10] == 6 };		// expected values
    for(int i = 0; i < (key_num-1)*2; i++){
      std::cout <<"nodes[" <<i<<"]"<<brt.internal_nodes[i]<<std::endl;
      std::cout <<"bitset: "<<   std::bitset<32>(brt.internal_nodes[i])<<std::endl;
    }
    binary_radix_tree::create_octree_nodes(key_num, brt);
		std::cout << "paper test " << ((allright) ? "ok." : "not ot.") << std::endl;
    };
int main() {


  Octomap::Octree<Point> octree = Octomap::Octree<Point>();

  std::vector<Point> pointcloud{};
  std::vector<uint32_t> sortedMorton{};
  std::vector<Code> mortonCodes{};
  for (int i = 0; i < 50; ++i) {
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
    mortonCodes.push_back(Code(mortonCode, 0));
    sortedMorton.push_back(mortonCode);
  }

    test_basic(sortedMorton);
  for(const Code& code: mortonCodes){
      octree.insertNode(code);
  }



  
/*
  
  auto end = std::chrono::high_resolution_clock::now();
  auto duration =std::chrono::duration_cast<std::chrono::milliseconds>(end - start)    .count();

  std::cout << "Time taken for insertion: " << duration << " milliseconds"     << std::endl;

  start = std::chrono::high_resolution_clock::now();
  for (auto point : pointcloud) {
    Code code = Code(encoder.morton3D(point), 0);
     std::cout <<"search code: "<< std::bitset<32>(code.getCode())<<std::endl;
    const Point *result = octree.search(code);
  }

  end = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start)              .count();

  std::cout << "Time taken for search: " << duration << " milliseconds"         << std::endl;
*/
  std::cout << "numer of leaf nodes: " << octree.getNumLeafNodes() << std::endl;
  std::cout << "numer of inner nodes: " << octree.getNumInnerNodes()
            << std::endl;
  std::cout << "numer of inner leaf nodes: " << octree.getNumInnerLeafNodes()
            << std::endl;
}