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
#include <execution>

#include "binary_radix_tree.hpp"
#include "morton_util.hpp"
#include "octree.hpp"
#include "util.hpp"

const int num_threads = std::thread::hardware_concurrency();
std::mutex mtx; // Mutex for protecting shared data

template <uint8_t Axis>
bool CompareAxis(const Eigen::Vector3f& a, const Eigen::Vector3f& b) {
  if constexpr (Axis == 0) {
    return a.x() < b.x();
  } else if constexpr (Axis == 1) {
    return a.y() < b.y();
  } else {
    return a.z() < b.z();
  }
}
// Function to perform counting sort on a specific digit's place (0 or 1)
void countingSort(std::vector<Code_t>& arr, int exp, int threadId, int numThreads) {
    const size_t n = arr.size();
    std::vector<Code_t> output(n);
    std::vector<Code_t> count(2, 0);

    for (size_t i = threadId; i < n; i += numThreads) {
        count[(arr[i] >> exp) & 1]++;
    }

    // Synchronize threads before updating the count array
    std::mutex mtx;
    std::unique_lock<std::mutex> lock(mtx);
    for (int i = 0; i < numThreads; i++) {
        count[0] += count[2 * i];
        count[1] += count[2 * i + 1];
    }
    lock.unlock();

    for (size_t i = threadId; i < n; i += numThreads) {
        output[count[(arr[i] >> exp) & 1]++] = arr[i];
    }

    for (size_t i = threadId; i < n; i += numThreads) {
        arr[i] = output[i];
    }
}

// Function to perform parallel radix sort using multiple threads
void parallelRadixSort(std::vector<Code_t>& arr, int numThreads) {
    //const size_t n = arr.size();
    const int numBits = 64; // Assuming 64-bit Morton codes

    for (int exp = 0; exp < numBits; exp++) {
        std::vector<std::thread> threads(numThreads);

        for (int i = 0; i < numThreads; i++) {
            threads[i] = std::thread(countingSort, std::ref(arr), exp, i, numThreads);
        }

        for (int i = 0; i < numThreads; i++) {
            threads[i].join();
        }
    }
}

void compute_morton_code_threaded(const int input_size, std::vector<Eigen::Vector3f>& inputs,std::vector<Code_t>& morton_keys, float min_coord, const float range ){

	const auto elements_per_thread = math::divide_ceil<int>(input_size, num_threads);
  std::cout<<"elements_per_thread: "<<elements_per_thread<<std::endl;
  	const auto worker_fn = [&inputs, &morton_keys,input_size, min_coord, range, elements_per_thread](int i) {
      int start = i * elements_per_thread;
      int end = math::min(input_size, (i + 1)*elements_per_thread);
      auto insertPosition = morton_keys.begin() + start;
      std::insert_iterator<std::vector<Code_t>> insertIter(morton_keys, insertPosition);
			          std::transform(inputs.begin() + start, inputs.begin() + end,
                            insertIter, [&](const auto& vec) {
                               return PointToCode(vec.x(), vec.y(), vec.z(), min_coord, range);
                           });
	};


	// Create the threads
	std::vector<std::thread> workers;
	for (int i = 0; i < num_threads; ++i)
		workers.push_back(std::thread(worker_fn, i));
	for (auto& t : workers)	t.join();
}

int main(int argc, char** argv) {
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

  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    return EXIT_SUCCESS;
  }

  if (!result.count("file")) {
    std::cerr
        << "requires an input file (\"../../data/1m_nn_uniform_4f.dat\")\n";
    std::cout << options.help() << std::endl;
    return EXIT_FAILURE;
  }

  const auto data_file = result["file"].as<std::string>();
  //const auto num_threads = result["thread"].as<int>();
  //const auto cpu = result["cpu"].as<bool>();
  const auto print = result["print"].as<bool>();

  thread_local std::mt19937 gen(114514);  // NOLINT(cert-msc51-cpp)
  static std::uniform_real_distribution dis(0.0f, 1.0f);

  // Prepare Inputs
  // constexpr int input_size = 1024;
  constexpr int input_size = 1280 * 720;
  std::vector<Eigen::Vector3f> inputs(input_size);
  std::generate(inputs.begin(), inputs.end(), [&] {
    const auto x = dis(gen) * 1024.0f;
    const auto y = dis(gen) * 1024.0f;
    const auto z = dis(gen) * 1024.0f;
    return Eigen::Vector3f(x, y, z);
  });

  float min_coord = 0.0f;
  float max_coord = 1.0f;

  TimeTask("Find Min Max", [&] {
    const auto [x_min, x_max] =
        std::minmax_element(inputs.begin(), inputs.end(), CompareAxis<0>);
    const auto [y_min, y_max] =
        std::minmax_element(inputs.begin(), inputs.end(), CompareAxis<1>);
    const auto [z_min, z_max] =
        std::minmax_element(inputs.begin(), inputs.end(), CompareAxis<2>);
    std::array<float, 3> mins{x_min->x(), y_min->y(), z_min->z()};
    std::array<float, 3> maxes{x_max->x(), y_max->y(), z_max->z()};
    min_coord = *std::min_element(mins.begin(), mins.end());
    max_coord = *std::max_element(maxes.begin(), maxes.end());
  });

  const float range = max_coord - min_coord;

  std::cout << "Min: " << min_coord << "\n";
  std::cout << "Max: " << max_coord << "\n";
  std::cout << "Range: " << range << "\n";

  // [Step 1] Compute Morton Codes
  std::vector<Code_t> morton_keys;
  morton_keys.reserve(input_size);

  
  TimeTask("Compute Morton Codes", [&] {
    std::transform(inputs.begin(), inputs.end(),
                   std::back_inserter(morton_keys), [&](const auto& vec) {
                     return PointToCode(vec.x(), vec.y(), vec.z(), min_coord,
                                        range);
                   });
  });
  
  //TimeTask("Compute Morton Codes", [&]{compute_morton_code_threaded(input_size, inputs, morton_keys, min_coord, range); });

  // [Step 2] Sort Morton Codes by Key
  TimeTask("Sort Morton Codes",
           [&] { std::sort(morton_keys.begin(), morton_keys.end()); });

  // [Step 3-4] Handle Duplicates
  TimeTask("Handle Duplicates", [&] {
    morton_keys.erase(std::unique(morton_keys.begin(), morton_keys.end()),
                      morton_keys.end());
  });

  if (print) {
    std::for_each(morton_keys.begin(), morton_keys.end(),
                  [min_coord, range](const auto key) {
                    std::cout << key << "\t" << std::bitset<kCodeLen>(key)
                              << "\t"
                              << CodeToPoint(key, min_coord, range).transpose()
                              << std::endl;
                  });
  }
  
  // [Step 5] Build Binary Radix Tree
  const auto num_brt_nodes = morton_keys.size() - 1;
  std::vector<brt::InnerNodes> inners(num_brt_nodes);

  TimeTask("Build Binary Radix Tree", [&] {
    create_binary_radix_tree_threaded(morton_keys.size(), morton_keys.data(), inners.data(), num_threads);
  });

  if (print) {
    for (unsigned int i = 0; i < num_brt_nodes; ++i) {
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
  TimeTask("Count Edges", [&] {
    // Copy a "1" to the first element to account for the root
    edge_count[0] = 1;
    oct::CalculateEdgeCount(edge_count.data(), inners.data(), num_brt_nodes);
  });

  // [Step 6.1] Prefix sum
  std::vector<int> oc_node_offsets(num_brt_nodes + 1);
  TimeTask("Prefix Sum", [&] {
    std::inclusive_scan(std::execution::par, edge_count.begin(), edge_count.end(),
                     oc_node_offsets.begin() + 1);
    oc_node_offsets[0] = 0;
  });

  // [Step 6.2] Allocate BH nodes
  const int num_oc_nodes = oc_node_offsets.back();
  const int root_level = inners[0].delta_node / 3;
  const Code_t root_prefix = morton_keys[0] >> (kCodeLen - (3 * root_level));
  std::vector<oct::OctNode> bh_nodes(num_oc_nodes);

  // Debug print
  std::cout << "Num Morton Keys: " << morton_keys.size() << "\n";
  std::cout << "Num Radix Nodes: " << num_brt_nodes << "\n";
  std::cout << "Num Octree Nodes: " << num_oc_nodes << "\n";

  // [Step 7] Create unlinked BH nodes
  TimeTask("Make Unlinked BH nodes", [&] {
    MakeNodesThreaded(bh_nodes.data(), oc_node_offsets.data(), edge_count.data(),
              morton_keys.data(), inners.data(), num_brt_nodes, num_threads, range);
  });

  // [Step 8] Linking BH nodes
  /*
  TimeTask("Link BH nodes", [&] {
    LinkNodes(bh_nodes.data(), oc_node_offsets.data(), edge_count.data(),
              morton_keys.data(), inners.data(), num_brt_nodes);
  });
  */
  
   TimeTask("Link BH nodes", [&] {
    LinkNodesThreaded(bh_nodes.data(), oc_node_offsets.data(), edge_count.data(),
              morton_keys.data(), inners.data(), num_brt_nodes, num_threads);
  });

  CheckTree(root_prefix, root_level * 3, bh_nodes.data(), 0,
            morton_keys.data());

  if (print) {
    for (int i = 0; i < num_oc_nodes; ++i) {
      std::cout << "OctNode " << i << "\n";
      std::cout << "\tchild_node_mask: "
                << std::bitset<8>(bh_nodes[i].child_node_mask) << "\n";
      std::cout << "\tchild_leaf_mask: "
                << std::bitset<8>(bh_nodes[i].child_leaf_mask) << "\n";

      std::cout << "\tchild : [" << bh_nodes[i].children[0];
      for (int j = 1; j < 8; ++j) {
        std::cout << ", " << bh_nodes[i].children[j];
      }
      std::cout << "]\n";

      std::cout << "\tcell_size: " << bh_nodes[i].cell_size << "\n";
      std::cout << "\tcornor: (" << bh_nodes[i].cornor.transpose() << ")\n";
      std::cout << "\n";
    }
  }


  return EXIT_SUCCESS;
}

