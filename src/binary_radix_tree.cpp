#include "include/binary_radix_tree.hpp"

#include <cassert>
#include <thread>

#include "include/morton_util.hpp"

namespace brt {

void process_internal_node(const int key_num, const Code_t* morton_keys,
                          InnerNodes* brt_nodes, int i) {
    const auto direction{
        math::sign<int>(Delta(morton_keys, i, i + 1) -
                        DeltaSafe(key_num, morton_keys, i, i - 1))};
    assert(direction != 0);

    const auto delta_min{DeltaSafe(key_num, morton_keys, i, i - direction)};

    int I_max{2};
    while (DeltaSafe(key_num, morton_keys, i, i + I_max * direction) >
           delta_min)
      I_max <<= 2;  // aka, *= 2

    // Find the other end using binary search.
    int I{0};
    for (int t{I_max / 2}; t; t /= 2)
      if (DeltaSafe(key_num, morton_keys, i, i + (I + t) * direction) >
          delta_min)
        I += t;

    const int j{i + I * direction};

    // Find the split position using binary search.
    const auto delta_node{DeltaSafe(key_num, morton_keys, i, j)};
    auto s{0};

    int t{I};
    do {
      t = math::divide2_ceil<int>(t);
      if (DeltaSafe(key_num, morton_keys, i, i + (s + t) * direction) >
          delta_node)
        s += t;
    } while (t > 1);

    const auto split{i + s * direction + math::min<int>(direction, 0)};

    // sanity check
    assert(Delta(morton_keys, i, j) > delta_min);
    assert(Delta(morton_keys, split, split + 1) == Delta(morton_keys, i, j));
    assert(!(split < 0 || split + 1 >= key_num));

    const int left{math::min<int>(i, j) == split
                       ? node::make_leaf<int>(split)
                       : node::make_internal<int>(split)};

    const int right{math::max<int>(i, j) == split + 1
                        ? node::make_leaf<int>(split + 1)
                        : node::make_internal<int>(split + 1)};

    brt_nodes[i].delta_node = delta_node;
    brt_nodes[i].left = left;
    brt_nodes[i].right = right;

    if (math::min<int>(i, j) != split) brt_nodes[left].parent = i;
    if (math::max<int>(i, j) != split + 1) brt_nodes[right].parent = i;
  
}
void create_binary_radix_tree_threaded(const int key_num, const Code_t* morton_keys,
                          InnerNodes* brt_nodes, int thread_number) {
	// compute the number of elements a thread will cover
	const auto elements_per_thread = math::divide_ceil<int>(key_num, thread_number);
	
	const auto worker_fn = [key_num, &morton_keys, &brt_nodes, elements_per_thread](int i) {
		for (int t = i * elements_per_thread; t < math::min(key_num - 1, (i + 1)*elements_per_thread); ++t)
			process_internal_node(key_num, morton_keys, brt_nodes, t);
	};

	// Create the threads
	std::vector<std::thread> workers;
	for (int i = 0; i < thread_number; ++i)
		workers.push_back(std::thread(worker_fn, i));
	for (auto& t : workers)	t.join();
}

void create_binary_radix_tree(const int key_num, const Code_t* morton_keys,
                          InnerNodes* brt_nodes) {
	for (int i = 0; i < key_num-1; ++i)
		process_internal_node(key_num, morton_keys, brt_nodes, i);
  }
}  // namespace brt