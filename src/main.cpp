#include <Eigen/Dense>
#include <thread>
#include <algorithm>
#include <array>
#include <bitset>
#include <cstdlib>
#include <cxxopts.hpp>
#include <iostream>
#include <numeric>
#include <mutex>
#include <random>
#include <omp.h>
#include <execution>
//#include <librealsense2/rs.hpp>
//#include <librealsense2/h/rs_types.h>
#include "binary_radix_tree.hpp"
#include "morton_util.hpp"
#include "octree.hpp"
#include "util.hpp"
// #include "occupancy_map.hpp"
using pcl_ptr = pcl::PointCloud<pcl::PointXYZ>::Ptr;
// const int num_threads = std::thread::hardware_concurrency();
std::mutex mtx; // Mutex for protecting shared data
/*
void test()
{
  float compute_time;
  float sort_time;
  float duplicate_time;
  float brt_time;
  float edges_time;
  float prefix_time;
  float make_nodes_time;
  float link_nodes_time;
  float total_time;

  int num_threads = 6;
  // Declare pointcloud object, for calculating pointclouds and texture mappings
  rs2::pointcloud pc;
  // We want the points object to be persistent so we can display the last cloud when a frame drops
  rs2::points points;

  // Declare RealSense pipeline, encapsulating the actual device and sensors
  rs2::pipeline pipe;
  // Start streaming with default recommended configuration
  pipe.start();
  oct::OctNode<float> *nodePtr = nullptr;

  std::vector<Code_t> old_pc;
  oct::Octree octree = oct::Octree(0, nodePtr, 0, nullptr);
  while (true)
  {
    // Wait for the next set of frames from the camera
    auto frames = pipe.wait_for_frames();

    auto depth = frames.get_depth_frame();

    points = pc.calculate(depth);
    auto cloud = points_to_pcl(points);
    auto inputs = cloud->points;
    int input_size = inputs.size();
    std::cout << "input size: " << input_size << std::endl;
    for (int i = 0; i < 50; ++i)
    {
      std::cout << inputs[i] << std::endl;
    }
    if (old_pc.empty())
    {

      float min_coord = 0.0f;
      float max_coord = 1.0f;
      TimeTask("Find Min Max", [&]
               {
    const auto [x_min, x_max] =
        std::minmax_element(inputs.begin(), inputs.end(), CompareAxis_pcl<0>);
    const auto [y_min, y_max] =
        std::minmax_element(inputs.begin(), inputs.end(), CompareAxis_pcl<1>);
    const auto [z_min, z_max] =
        std::minmax_element(inputs.begin(), inputs.end(), CompareAxis_pcl<2>);
    std::array<float, 3> mins{x_min->x, y_min->y, z_min->z};
    std::array<float, 3> maxes{x_max->x, y_max->y, z_max->z};
    min_coord = *std::min_element(mins.begin(), mins.end());
    max_coord = *std::max_element(maxes.begin(), maxes.end()); });

      const float range = max_coord - min_coord;

      std::cout << "Min: " << min_coord << "\n";
      std::cout << "Max: " << max_coord << "\n";
      std::cout << "Range: " << range << "\n";
      // [Step 1] Compute Morton Codes
      std::vector<Code_t> morton_keys(input_size, 0xffffffff);
      // morton_keys.reserve(input_size);

      compute_time = TimeTask("Compute Morton Codes", [&]
                              { compute_morton_code_openmp_pcl(input_size, inputs, morton_keys, min_coord, range, num_threads); });

      // [Step 2] Sort Morton Codes by Key
      sort_time = TimeTask("Sort Morton Codes",
                           [&]
                           { std::sort(std::execution::par, morton_keys.begin(), morton_keys.end()); });

      sort_time = TimeTask("Sort Morton Codes radix sort",
                           [&]
                           { omp_lsd_radix_sort(morton_keys.size(), morton_keys, num_threads); });
      // [Step 3-4] Handle Duplicates
      duplicate_time = TimeTask("Handle Duplicates", [&]
                                { morton_keys.erase(std::unique(morton_keys.begin(), morton_keys.end()),
                                                    morton_keys.end()); });

      // [Step 5] Build Binary Radix Tree
      const auto num_brt_nodes = morton_keys.size() - 1;
      std::vector<brt::InnerNodes> inners(num_brt_nodes);

      brt_time = TimeTask("Build Binary Radix Tree", [&]
                          { create_binary_radix_tree_threaded(morton_keys.size(), morton_keys.data(), inners.data(), num_threads); });

      // [Step 6] Count edges
      std::vector<int> edge_count(num_brt_nodes);

      edges_time = TimeTask("Count Edges", [&]
                            {
    // Copy a "1" to the first element to account for the root
    edge_count[0] = 1;
    oct::CalculateEdgeCountThreaded(edge_count.data(), inners.data(), num_brt_nodes, num_threads); });

      // [Step 6.1] Prefix sum

      std::vector<int> oc_node_offsets(num_brt_nodes + 1);

      prefix_time = TimeTask("Prefix Sum", [&]
                             {
    PrefixSumThreaded(edge_count, oc_node_offsets, num_brt_nodes, num_threads);
    oc_node_offsets[0] = 0; });

      // [Step 6.2] Allocate BH nodes
      const int num_oc_nodes = oc_node_offsets.back();
      const int root_level = inners[0].delta_node / 3;
      const Code_t root_prefix = morton_keys[0] >> (kCodeLen - (3 * root_level));
      std::vector<oct::OctNode<float>> bh_nodes(num_oc_nodes);

      // Debug print
      std::cout << "Num Morton Keys: " << morton_keys.size() << "\n";
      std::cout << "Num Radix Nodes: " << num_brt_nodes << "\n";
      std::cout << "Num Octree Nodes: " << num_oc_nodes << "\n";

      // [Step 7] Create unlinked BH nodes
      make_nodes_time = TimeTask("Make Unlinked BH nodes", [&]
                                 { oct::MakeNodesThreaded(bh_nodes.data(), oc_node_offsets.data(), edge_count.data(),
                                                          morton_keys.data(), inners.data(), num_brt_nodes, num_threads, range); });

      // [Step 8] Linking BH nodes
      link_nodes_time = TimeTask("Link BH nodes", [&]
                                 { oct::LinkNodesThreaded(bh_nodes.data(), oc_node_offsets.data(), edge_count.data(),
                                                          morton_keys.data(), inners.data(), num_brt_nodes, num_threads); });

      octree = oct::Octree(root_prefix, bh_nodes.data(), root_level, morton_keys.data());
      old_pc = morton_keys;
    }
    else
    {
      float min_coord = 0.0f;
      float max_coord = 1.0f;
      TimeTask("Find Min Max", [&]
               {
    const auto [x_min, x_max] =
        std::minmax_element(inputs.begin(), inputs.end(), CompareAxis_pcl<0>);
    const auto [y_min, y_max] =
        std::minmax_element(inputs.begin(), inputs.end(), CompareAxis_pcl<1>);
    const auto [z_min, z_max] =
        std::minmax_element(inputs.begin(), inputs.end(), CompareAxis_pcl<2>);
    std::array<float, 3> mins{x_min->x, y_min->y, z_min->z};
    std::array<float, 3> maxes{x_max->x, y_max->y, z_max->z};
    min_coord = *std::min_element(mins.begin(), mins.end());
    max_coord = *std::max_element(maxes.begin(), maxes.end()); });

      const float range = max_coord - min_coord;
      std::vector<Code_t> morton_keys(input_size, 0xffffffff);
      compute_morton_code_openmp_pcl(input_size, inputs, morton_keys, min_coord, range, num_threads);

      std::vector<Code_t> to_add = octree.update_nodes(morton_keys);
      if (!to_add.empty())
      {
        morton_keys.insert(morton_keys.end(), to_add.begin(), to_add.end());
        std::cout << "reconstructing" << std::endl;
        // [Step 2] Sort Morton Codes by Key
        sort_time = TimeTask("Sort Morton Codes",
                             [&]
                             { std::sort(std::execution::par, morton_keys.begin(), morton_keys.end()); });

        sort_time = TimeTask("Sort Morton Codes radix sort",
                             [&]
                             { omp_lsd_radix_sort(morton_keys.size(), morton_keys, num_threads); });
        // [Step 3-4] Handle Duplicates
        duplicate_time = TimeTask("Handle Duplicates", [&]
                                  { morton_keys.erase(std::unique(morton_keys.begin(), morton_keys.end()),
                                                      morton_keys.end()); });

        // [Step 5] Build Binary Radix Tree
        const auto num_brt_nodes = morton_keys.size() - 1;
        std::vector<brt::InnerNodes> inners(num_brt_nodes);

        brt_time = TimeTask("Build Binary Radix Tree", [&]
                            { create_binary_radix_tree_threaded(morton_keys.size(), morton_keys.data(), inners.data(), num_threads); });

        // [Step 6] Count edges
        std::vector<int> edge_count(num_brt_nodes);

        edges_time = TimeTask("Count Edges", [&]
                              {
    // Copy a "1" to the first element to account for the root
    edge_count[0] = 1;
    oct::CalculateEdgeCountThreaded(edge_count.data(), inners.data(), num_brt_nodes, num_threads); });

        // [Step 6.1] Prefix sum

        std::vector<int> oc_node_offsets(num_brt_nodes + 1);

        prefix_time = TimeTask("Prefix Sum", [&]
                               {
    PrefixSumThreaded(edge_count, oc_node_offsets, num_brt_nodes, num_threads);
    oc_node_offsets[0] = 0; });

        // [Step 6.2] Allocate BH nodes
        const int num_oc_nodes = oc_node_offsets.back();
        const int root_level = inners[0].delta_node / 3;
        const Code_t root_prefix = morton_keys[0] >> (kCodeLen - (3 * root_level));
        std::vector<oct::OctNode<float>> bh_nodes(num_oc_nodes);

        // Debug print
        std::cout << "Num Morton Keys: " << morton_keys.size() << "\n";
        std::cout << "Num Radix Nodes: " << num_brt_nodes << "\n";
        std::cout << "Num Octree Nodes: " << num_oc_nodes << "\n";

        // [Step 7] Create unlinked BH nodes
        make_nodes_time = TimeTask("Make Unlinked BH nodes", [&]
                                   { oct::MakeNodesThreaded(bh_nodes.data(), oc_node_offsets.data(), edge_count.data(),
                                                            morton_keys.data(), inners.data(), num_brt_nodes, num_threads, range); });

        // [Step 8] Linking BH nodes
        link_nodes_time = TimeTask("Link BH nodes", [&]
                                   { oct::LinkNodesThreaded(bh_nodes.data(), oc_node_offsets.data(), edge_count.data(),
                                                            morton_keys.data(), inners.data(), num_brt_nodes, num_threads); });

        octree = oct::Octree(root_prefix, bh_nodes.data(), root_level, morton_keys.data());
        old_pc = morton_keys;
      }
    }
  }
}


*/


void test2()
{
  float compute_time;
  float sort_time;
  float duplicate_time;
  float brt_time;
  float edges_time;
  float prefix_time;
  float make_nodes_time;
  float link_nodes_time;
  float total_time;

  int num_threads = 6;
  //oct::OctNode<float> *nodePtr = nullptr;
  std::vector<oct::OctNode<float>> emptyoctnode;
  std::vector<Code_t> emptycode;
  std::vector<Code_t> old_pc;
  oct::Octree octree = oct::Octree(0, emptyoctnode, 0, emptycode);
    thread_local std::mt19937 gen(114514); // NOLINT(cert-msc51-cpp)
  static std::uniform_real_distribution dis(0.0f, 1.0f);

  // Prepare Inputs
  // constexpr int input_size = 1024;
  constexpr int input_size = 1280 * 720;
  std::vector<Eigen::Vector3f> inputs(input_size);
  for(int j = 0; j < 10; ++j)
  {
   std::generate(inputs.begin(), inputs.end(), [&]
                {
    const auto x = dis(gen) * 1024.0f;
    const auto y = dis(gen) * 1024.0f;
    const auto z = dis(gen) * 1024.0f;
    return Eigen::Vector3f(x, y, z); });
    std::cout << "input size: " << input_size << std::endl;

    if (old_pc.empty())
    {
      std::cout<<"empty, construct octree first"<<std::endl;

      float min_coord = 0.0f;
      float max_coord = 1.0f;
      TimeTask("Find Min Max", [&]
               {
    const auto [x_min, x_max] =
        std::minmax_element(inputs.begin(), inputs.end(), CompareAxis<0>);
    const auto [y_min, y_max] =
        std::minmax_element(inputs.begin(), inputs.end(), CompareAxis<1>);
    const auto [z_min, z_max] =
        std::minmax_element(inputs.begin(), inputs.end(), CompareAxis<2>);
    std::array<float, 3> mins{x_min->x(), y_min->y(), z_min->z()};
    std::array<float, 3> maxes{x_max->x(), y_max->y(), z_max->z()};
    min_coord = *std::min_element(mins.begin(), mins.end());
    max_coord = *std::max_element(maxes.begin(), maxes.end()); });

      const float range = max_coord - min_coord;

      std::cout << "Min: " << min_coord << "\n";
      std::cout << "Max: " << max_coord << "\n";
      std::cout << "Range: " << range << "\n";
      // [Step 1] Compute Morton Codes
      std::vector<Code_t> morton_keys(input_size, 0xffffffff);
      // morton_keys.reserve(input_size);

      compute_time = TimeTask("Compute Morton Codes", [&]
                              { compute_morton_code_openmp(input_size, inputs, morton_keys, min_coord, range, num_threads); });

      // [Step 2] Sort Morton Codes by Key
      sort_time = TimeTask("Sort Morton Codes",
                           [&]
                           { std::sort(std::execution::par, morton_keys.begin(), morton_keys.end()); });

      // [Step 3-4] Handle Duplicates
      duplicate_time = TimeTask("Handle Duplicates", [&]
                                { morton_keys.erase(std::unique(morton_keys.begin(), morton_keys.end()),
                                                    morton_keys.end()); });

      // [Step 5] Build Binary Radix Tree
      const auto num_brt_nodes = morton_keys.size() - 1;
      std::vector<brt::InnerNodes> inners(num_brt_nodes);

      brt_time = TimeTask("Build Binary Radix Tree", [&]
                          { create_binary_radix_tree_threaded(morton_keys.size(), morton_keys.data(), inners.data(), num_threads); });

      // [Step 6] Count edges
      std::vector<int> edge_count(num_brt_nodes);

      edges_time = TimeTask("Count Edges", [&]
                            {
    // Copy a "1" to the first element to account for the root
    edge_count[0] = 1;
    oct::CalculateEdgeCountThreaded(edge_count.data(), inners.data(), num_brt_nodes, num_threads); });

      // [Step 6.1] Prefix sum

      std::vector<int> oc_node_offsets(num_brt_nodes + 1);

      prefix_time = TimeTask("Prefix Sum", [&]
                             {
    PrefixSumThreaded(edge_count, oc_node_offsets, num_brt_nodes, num_threads);
    oc_node_offsets[0] = 0; });

      // [Step 6.2] Allocate BH nodes
      const int num_oc_nodes = oc_node_offsets.back();
      const int root_level = inners[0].delta_node / 3;
      const Code_t root_prefix = morton_keys[0] >> (kCodeLen - (3 * root_level));
      std::vector<oct::OctNode<float>> bh_nodes(num_oc_nodes);

      // Debug print
      std::cout << "Num Morton Keys: " << morton_keys.size() << "\n";
      std::cout << "Num Radix Nodes: " << num_brt_nodes << "\n";
      std::cout << "Num Octree Nodes: " << num_oc_nodes << "\n";

      // [Step 7] Create unlinked BH nodes
      make_nodes_time = TimeTask("Make Unlinked BH nodes", [&]
                                 { oct::MakeNodesThreaded(bh_nodes.data(), oc_node_offsets.data(), edge_count.data(),
                                                          morton_keys.data(), inners.data(), num_brt_nodes, num_threads, range); });

      // [Step 8] Linking BH nodes
      link_nodes_time = TimeTask("Link BH nodes", [&]
                                 { oct::LinkNodesThreaded(bh_nodes.data(), oc_node_offsets.data(), edge_count.data(),
                                                          morton_keys.data(), inners.data(), num_brt_nodes, num_threads); });

      octree = oct::Octree(root_prefix, bh_nodes, root_level, morton_keys);
      old_pc = morton_keys;
    }
    else
    {
      float min_coord = 0.0f;
      float max_coord = 1.0f;
      TimeTask("Find Min Max", [&]
               {
    const auto [x_min, x_max] =
        std::minmax_element(inputs.begin(), inputs.end(), CompareAxis<0>);
    const auto [y_min, y_max] =
        std::minmax_element(inputs.begin(), inputs.end(), CompareAxis<1>);
    const auto [z_min, z_max] =
        std::minmax_element(inputs.begin(), inputs.end(), CompareAxis<2>);
    std::array<float, 3> mins{x_min->x(), y_min->y(), z_min->z()};
    std::array<float, 3> maxes{x_max->x(), y_max->y(), z_max->z()};
    min_coord = *std::min_element(mins.begin(), mins.end());
    max_coord = *std::max_element(maxes.begin(), maxes.end()); });

      const float range = max_coord - min_coord;
      std::vector<Code_t> morton_keys(input_size, 0xffffffff);
      compute_morton_code_openmp(input_size, inputs, morton_keys, min_coord, range, num_threads);

      std::vector<Code_t> to_add;
      TimeTask("Update Node", [&]{to_add = octree.update_nodes(morton_keys);});
      if (!to_add.empty())
      {
        morton_keys.insert(morton_keys.end(), to_add.begin(), to_add.end());
        std::cout << "reconstructing" << std::endl;
        // [Step 2] Sort Morton Codes by Key
        sort_time = TimeTask("Sort Morton Codes",
                             [&]
                             { std::sort(std::execution::par, morton_keys.begin(), morton_keys.end()); });
        // [Step 3-4] Handle Duplicates
        duplicate_time = TimeTask("Handle Duplicates", [&]
                                  { morton_keys.erase(std::unique(morton_keys.begin(), morton_keys.end()),
                                                      morton_keys.end()); });

        // [Step 5] Build Binary Radix Tree
        const auto num_brt_nodes = morton_keys.size() - 1;
        std::vector<brt::InnerNodes> inners(num_brt_nodes);

        brt_time = TimeTask("Build Binary Radix Tree", [&]
                            { create_binary_radix_tree_threaded(morton_keys.size(), morton_keys.data(), inners.data(), num_threads); });

        // [Step 6] Count edges
        std::vector<int> edge_count(num_brt_nodes);

        edges_time = TimeTask("Count Edges", [&]
                              {
    // Copy a "1" to the first element to account for the root
    edge_count[0] = 1;
    oct::CalculateEdgeCountThreaded(edge_count.data(), inners.data(), num_brt_nodes, num_threads); });

        // [Step 6.1] Prefix sum

        std::vector<int> oc_node_offsets(num_brt_nodes + 1);

        prefix_time = TimeTask("Prefix Sum", [&]
                               {
    PrefixSumThreaded(edge_count, oc_node_offsets, num_brt_nodes, num_threads);
    oc_node_offsets[0] = 0; });

        // [Step 6.2] Allocate BH nodes
        const int num_oc_nodes = oc_node_offsets.back();
        const int root_level = inners[0].delta_node / 3;
        const Code_t root_prefix = morton_keys[0] >> (kCodeLen - (3 * root_level));
        std::vector<oct::OctNode<float>> bh_nodes(num_oc_nodes);

        // Debug print
        std::cout << "Num Morton Keys: " << morton_keys.size() << "\n";
        std::cout << "Num Radix Nodes: " << num_brt_nodes << "\n";
        std::cout << "Num Octree Nodes: " << num_oc_nodes << "\n";

        // [Step 7] Create unlinked BH nodes
        make_nodes_time = TimeTask("Make Unlinked BH nodes", [&]
                                   { oct::MakeNodesThreaded(bh_nodes.data(), oc_node_offsets.data(), edge_count.data(),
                                                            morton_keys.data(), inners.data(), num_brt_nodes, num_threads, range); });

        // [Step 8] Linking BH nodes
        link_nodes_time = TimeTask("Link BH nodes", [&]
                                   { oct::LinkNodesThreaded(bh_nodes.data(), oc_node_offsets.data(), edge_count.data(),
                                                            morton_keys.data(), inners.data(), num_brt_nodes, num_threads); });

        octree = oct::Octree(root_prefix, bh_nodes, root_level, morton_keys);
        old_pc = morton_keys;
      }
    }
  }
}

int main(int argc, char **argv)
{
  test2();
  float compute_time;
  float sort_time;
  float duplicate_time;
  float brt_time;
  float edges_time;
  float prefix_time;
  float make_nodes_time;
  float link_nodes_time;
  float total_time;
  cxxopts::Options options(
      "Redwood Radix Tree",
      "Redwood accelerated Binary Radix Tree and Octree Implementation");

  // clang-format off
  options.add_options()
    ("f,file", "Input file name", cxxopts::value<std::string>())
    ("j,thread", "Number of threads", cxxopts::value<int>()->default_value("1"))
    ("c,cpu", "Enable CPU baseline", cxxopts::value<bool>()->default_value("false"))
    ("p,print", "Print Debug", cxxopts::value<bool>()->default_value("false"))
    ("h,help", "Print usage");
  // clang-format on

  options.parse_positional({"file"});

  const auto result = options.parse(argc, argv);

  if (result.count("help"))
  {
    std::cout << options.help() << std::endl;
    return EXIT_SUCCESS;
  }

  if (!result.count("file"))
  {
    std::cerr
        << "requires an input file (\"../../data/1m_nn_uniform_4f.dat\")\n";
    std::cout << options.help() << std::endl;
    return EXIT_FAILURE;
  }

  const auto data_file = result["file"].as<std::string>();
  // const auto cpu = result["cpu"].as<bool>();
  const auto print = result["print"].as<bool>();
  const auto input_threads = result["thread"].as<int>();
  const int max_threads = std::thread::hardware_concurrency();
  const auto num_threads = std::min(input_threads, max_threads);
  std::cout << "num of threads: " << num_threads << std::endl;

  thread_local std::mt19937 gen(114514); // NOLINT(cert-msc51-cpp)
  static std::uniform_real_distribution dis(0.0f, 1.0f);

  // Prepare Inputs
  // constexpr int input_size = 1024;
  constexpr int input_size = 1280 * 720;
  std::vector<Eigen::Vector3f> inputs(input_size);
  std::generate(inputs.begin(), inputs.end(), [&]
                {
    const auto x = dis(gen) * 1024.0f;
    const auto y = dis(gen) * 1024.0f;
    const auto z = dis(gen) * 1024.0f;
    return Eigen::Vector3f(x, y, z); });
  float min_coord = 0.0f;
  float max_coord = 1.0f;

  TimeTask("Find Min Max", [&]
           {
    const auto [x_min, x_max] =
        std::minmax_element(inputs.begin(), inputs.end(), CompareAxis<0>);
    const auto [y_min, y_max] =
        std::minmax_element(inputs.begin(), inputs.end(), CompareAxis<1>);
    const auto [z_min, z_max] =
        std::minmax_element(inputs.begin(), inputs.end(), CompareAxis<2>);
    std::array<float, 3> mins{x_min->x(), y_min->y(), z_min->z()};
    std::array<float, 3> maxes{x_max->x(), y_max->y(), z_max->z()};
    min_coord = *std::min_element(mins.begin(), mins.end());
    max_coord = *std::max_element(maxes.begin(), maxes.end()); });

  const float range = max_coord - min_coord;

  std::cout << "Min: " << min_coord << "\n";
  std::cout << "Max: " << max_coord << "\n";
  std::cout << "Range: " << range << "\n";

  // [Step 1] Compute Morton Codes
  std::vector<Code_t> morton_keys(input_size, 0xffffffff);
  // morton_keys.reserve(input_size);

  compute_time = TimeTask("Compute Morton Codes", [&]
                          { compute_morton_code_openmp(input_size, inputs, morton_keys, min_coord, range, num_threads); });
  sort_time = TimeTask("Sort Morton Codes",
                       [&]
                       { std::sort(std::execution::par, morton_keys.begin(), morton_keys.end()); });

  // [Step 3-4] Handle Duplicates
  duplicate_time = TimeTask("Handle Duplicates", [&]
                            { morton_keys.erase(std::unique(morton_keys.begin(), morton_keys.end()),
                                                morton_keys.end()); });

  if (print)
  {
    std::for_each(morton_keys.begin(), morton_keys.end(),
                  [min_coord, range](const auto key)
                  {
                    std::cout << key << "\t" << std::bitset<kCodeLen>(key)
                              << "\t"
                              << CodeToPoint(key, min_coord, range).transpose()
                              << std::endl;
                  });
  }

  // [Step 5] Build Binary Radix Tree
  const auto num_brt_nodes = morton_keys.size() - 1;
  std::vector<brt::InnerNodes> inners(num_brt_nodes);

  brt_time = TimeTask("Build Binary Radix Tree", [&]
                      { create_binary_radix_tree_threaded(morton_keys.size(), morton_keys.data(), inners.data(), num_threads); });

  if (print)
  {
    for (unsigned int i = 0; i < num_brt_nodes; ++i)
    {
      std::cout << "Node " << i << "\n";
      std::cout << "\tdelta_node: " << inners[i].delta_node << "\n";
      std::cout << "\tleft: " << inners[i].left << "\n";
      std::cout << "\tright: " << inners[i].right << "\n";
      std::cout << "\tparent: " << inners[i].parent << "\n";
      std::cout << "\n";
    }
  }

  // [Step 6] Count edges
  std::vector<int> edge_count(num_brt_nodes);

  edges_time = TimeTask("Count Edges", [&]
                        {
    // Copy a "1" to the first element to account for the root
    edge_count[0] = 1;
    oct::CalculateEdgeCountThreaded(edge_count.data(), inners.data(), num_brt_nodes, num_threads); });

  // [Step 6.1] Prefix sum

  std::vector<int> oc_node_offsets(num_brt_nodes + 1);

  prefix_time = TimeTask("Prefix Sum", [&]
                         {
    PrefixSumThreaded(edge_count, oc_node_offsets, num_brt_nodes, num_threads);
    oc_node_offsets[0] = 0; });

  // [Step 6.2] Allocate BH nodes
  const int num_oc_nodes = oc_node_offsets.back();
  const int root_level = inners[0].delta_node / 3;
  const Code_t root_prefix = morton_keys[0] >> (kCodeLen - (3 * root_level));
  std::vector<oct::OctNode<float>> bh_nodes(num_oc_nodes);

  // Debug print
  std::cout << "Num Morton Keys: " << morton_keys.size() << "\n";
  std::cout << "Num Radix Nodes: " << num_brt_nodes << "\n";
  std::cout << "Num Octree Nodes: " << num_oc_nodes << "\n";

  // [Step 7] Create unlinked BH nodes
  make_nodes_time = TimeTask("Make Unlinked BH nodes", [&]
                             { oct::MakeNodesThreaded(bh_nodes.data(), oc_node_offsets.data(), edge_count.data(),
                                                      morton_keys.data(), inners.data(), num_brt_nodes, num_threads, range); });

  // [Step 8] Linking BH nodes
  link_nodes_time = TimeTask("Link BH nodes", [&]
                             { oct::LinkNodesThreaded(bh_nodes.data(), oc_node_offsets.data(), edge_count.data(),
                                                      morton_keys.data(), inners.data(), num_brt_nodes, num_threads); });

  oct::Octree octree = oct::Octree(root_prefix, bh_nodes, root_level, morton_keys);
  auto oct_idx = octree.search(morton_keys[800], root_prefix, root_level * 3, 0);
  std::cout << "oct index: " << oct_idx << std::endl;
  oct::CheckTree(root_prefix, root_level * 3, bh_nodes.data(), 0,
                 morton_keys.data());

  if (print)
  {
    for (int i = 0; i < num_oc_nodes; ++i)
    {
      std::cout << "OctNode " << i << "\n";
      std::cout << "\tchild_node_mask: "
                << std::bitset<8>(bh_nodes[i].child_node_mask) << "\n";
      std::cout << "\tchild_leaf_mask: "
                << std::bitset<8>(bh_nodes[i].child_leaf_mask) << "\n";

      std::cout << "\tchild : [" << bh_nodes[i].children[0];
      for (int j = 1; j < 8; ++j)
      {
        std::cout << ", " << bh_nodes[i].children[j];
      }
      std::cout << "]\n";

      std::cout << "\tcell_size: " << bh_nodes[i].cell_size << "\n";
      std::cout << "\tcornor: (" << bh_nodes[i].cornor.transpose() << ")\n";
      std::cout << "\n";
    }
  }

  total_time = compute_time +
               sort_time +
               duplicate_time +
               brt_time +
               edges_time +
               prefix_time +
               make_nodes_time +
               link_nodes_time;

  std::cout << "total time: " << total_time << std::endl;
  std::cout << "Compute Morton Codes Time: " << compute_time << " (" << compute_time / total_time * 100 << "%)" << std::endl;
  std::cout << "Sort Morton Codes Time: " << sort_time << " (" << sort_time / total_time * 100 << "%)" << std::endl;
  std::cout << "Handle Duplicates Time: " << duplicate_time << " (" << duplicate_time / total_time * 100 << "%)" << std::endl;
  std::cout << "Build Binary Radix Tree Time: " << brt_time << " (" << brt_time / total_time * 100 << "%)" << std::endl;
  std::cout << "Count Edges Time: " << edges_time << " (" << edges_time / total_time * 100 << "%)" << std::endl;
  std::cout << "Prefix Sum Time: " << prefix_time << " (" << prefix_time / total_time * 100 << "%)" << std::endl;
  std::cout << "Make Unlinked BH nodes Time: " << make_nodes_time << " (" << make_nodes_time / total_time * 100 << "%)" << std::endl;
  std::cout << "Link BH nodes Time: " << link_nodes_time << " (" << link_nodes_time / total_time * 100 << "%)" << std::endl;

  return EXIT_SUCCESS;
}
