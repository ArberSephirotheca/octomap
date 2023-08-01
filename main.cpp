#include <iostream>
#include "OctreeBase.cpp"
int main(){
    Octomap::Octree octree = Octomap::Octree();
    Point one = Point{0.2,0.4,0.5};
    Point two = Point{0.3, 0.5, 0.2};
    octree.insert(one);
    octree.insert(two);

    auto result1 = octree.query(one);
    for (const auto& data : result1) {
        std::cout << "Data: " << data << std::endl;
    }

    auto result2 = octree.query(two);
    for (const auto& data : result2) {
        std::cout << "Data: " << data << std::endl;
    }

    return 0;
}