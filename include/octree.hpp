#pragma once

#include <Eigen/Dense>

#include "binary_radix_tree.hpp"
#include <thread>
#include <bitset>
using Code_t = uint64_t;

namespace oct {

template<typename DATA_TYPE>
struct OctNode {
  // Payload
  DATA_TYPE data;

  Eigen::Vector3f cornor;
  float cell_size;

  // TODO: This is overkill number of pointers
  int children[8];

  /**
   * @brief For bit position i (from the right): If 1, children[i] is the index
   * of a child octree node. If 0, the ith child is either absent, or
   * children[i] is the index of a leaf.
   */
  int child_node_mask;

  /**
   * @brief For bit position i (from the right): If 1, children[i] is the index
   * of a leaf (in the corresponding points array). If 0, the ith child is
   * either absent, or an octree node.
   */
  int child_leaf_mask;

  /**
   * @brief Set a child
   *
   * @param child: index of octree node that will become the child
   * @param my_child_idx: which of my children it will be [0-7]
   */
  void SetChild(const int child, const int my_child_idx);

  /**
   * @brief Set the Leaf object
   *
   * @param leaf: index of point that will become the leaf child
   * @param my_child_idx: which of my children it will be [0-7]
   */
  void SetLeaf(const int leaf, const int my_child_idx);
};

template <typename DATA_TYPE>
class Octree{
  public:
  Octree(Code_t root_prefix, std::vector<OctNode<DATA_TYPE>> nodes, int level, std::vector<Code_t> codes): 
        root_prefix_(root_prefix), nodes_(nodes), level_(level), codes_(codes) {}
  ~Octree() {}
  Octree& operator=(const Octree& other){
    root_prefix_ = other.root_prefix_;
    nodes_ = other.nodes_;
    level_ = other.level_;
    codes_ = other.codes_;
    return *this;
  }
  int search(const Code_t code, const Code_t prefix, const int code_len, const int oct_idx);
  void update_node(const Code_t code);
  std::vector<Code_t> update_nodes(const std::vector<Code_t> codes);
  private:
    Code_t root_prefix_;
    std::vector<OctNode<DATA_TYPE>> nodes_;
    int level_;
    std::vector<Code_t> codes_;
    std::vector<Code_t> to_add;
};


template<typename DATA_TYPE>
void OctNode<DATA_TYPE>::SetChild(const int child, const int my_child_idx) {
  children[my_child_idx] = child;
  // TODO: atomicOr in CUDA
  child_node_mask |= (1 << my_child_idx);
}

template<typename DATA_TYPE>
void OctNode<DATA_TYPE>::SetLeaf(const int leaf, const int my_child_idx) {
  children[my_child_idx] = leaf;
  // TODO: atomicOr in CUDA
  child_leaf_mask |= (1 << my_child_idx);
}

template<typename DATA_TYPE>
  int Octree<DATA_TYPE>::search(const Code_t code, const Code_t prefix, int code_len, int oct_idx) {
    const int k_shift = kCodeLen - 3;
    size_t child_idx;
  while(code_len <= 63){
    OctNode<DATA_TYPE>& node = nodes_[oct_idx];
    child_idx = code >> (k_shift - code_len) &0x7;
  if(node.child_node_mask & (1 << child_idx)){
    code_len += 3;
    oct_idx = node.children[child_idx];
  } else if (node.child_leaf_mask & (1 << child_idx)){
    Code_t leaf_code = codes_[node.children[child_idx]];
    if(code == leaf_code){
      return node.children[child_idx];
    } else {
      return -1;
    }
  } else {
    return -1;
  }
  }
}

template<typename DATA_TYPE>
  void Octree<DATA_TYPE>::update_node(const Code_t code){
    volatile int oct_idx = search(code, root_prefix_, level_*3, 0);
    if(oct_idx == -1){
      to_add.push_back(code);
    }else{
      nodes_[oct_idx].data += 1;
    }
    
  }

template<typename DATA_TYPE>
  std::vector<Code_t> Octree<DATA_TYPE>::update_nodes(const std::vector<Code_t> codes){
    //parallize the lopp with openmp
    //#pragma omp parallel for
    for(int i = 0; i < codes.size(); i++){
      update_node(codes[i]);
    }
    return to_add;
  }


 bool IsLeaf(const int internal_value) {
  // check the most significant bit, which is used as a flag for "is leaf node"
  return internal_value >> (sizeof(int) * 8 - 1);
}

 int GetLeafIndex(const int internal_value) {
  // delete the last bit which tells if this is leaf or internal index
  return internal_value &
         ~(1 << (sizeof(int) * 8 -
                 1));  // NOLINT(clang-diagnostic-shift-sign-overflow)
}

/**
 * @brief Count the number of Octree node under a binary radix tree node.
 *
 * @param edge_count: results edge count array
 * @param inners: binary radix tree nodes
 * @param num_brt_nodes: number of binary radix tree nodes
 */

void CalculateEdgeCountThreaded(int* edge_count, const brt::InnerNodes* inners,
                        const int num_brt_nodes, int num_threads){
	const auto elements_per_thread = math::divide_ceil<int>(num_brt_nodes, num_threads);
	const auto worker_fn = [edge_count, inners, num_brt_nodes, elements_per_thread](int i) {
		for (int t = (i * elements_per_thread)+1; t < math::min(num_brt_nodes, (i + 1)*elements_per_thread+1); ++t){
        const int my_depth = inners[t].delta_node / 3;
        const int parent_depth = inners[inners[t].parent].delta_node / 3;
        edge_count[t] = my_depth - parent_depth;
      }
	};

	// Create the threads
	std::vector<std::thread> workers;
	for (int i = 0; i < num_threads; ++i)
		workers.push_back(std::thread(worker_fn, i));
	for (auto& t : workers)	t.join();

}
/**
 * @brief Make the unlinked octree nodes from the binary radix tree.
 * https://github.com/ahmidou/ShapeExtraction/blob/master/src/Octree.cu
 *
 * @param nodes: array of preallocated octree nodes
 * @param node_offsets: ranges of each RT node
 * @param edge_count: number of nodes in each RT node
 * @param morton_keys: sorted morton keys
 * @param inners: binary radix tree nodes
 * @param num_brt_nodes: number of binary radix tree nodes
 * @param tree_range: range of the entire octree, default 1.0f
 */
template<typename DATA_TYPE>
void MakeNodes(OctNode<DATA_TYPE>* nodes, const int* node_offsets, const int* edge_count,
               const Code_t* morton_keys, const brt::InnerNodes* inners,
               int num_brt_nodes, float tree_range = 1.0f);

template<typename DATA_TYPE>
void MakeNodesThreaded(OctNode<DATA_TYPE>* nodes, const int* node_offsets, const int* edge_count,
               const Code_t* morton_keys, const brt::InnerNodes* inners,
               int num_brt_nodes, int num_threads, float tree_range = 1.0f){
                
  // the root doesn't represent level 0 of the "entire" octree
  const int root_level = inners[0].delta_node / 3;
  const Code_t root_prefix = morton_keys[0] >> (kCodeLen - (root_level * 3));

  nodes[0].cornor = CodeToPoint(root_prefix << (kCodeLen - (root_level * 3)));
  nodes[0].cell_size = tree_range;

	const auto elements_per_thread = math::divide_ceil<int>(num_brt_nodes, num_threads);
  // skipping root

  	const auto worker_fn = [nodes, node_offsets, edge_count,morton_keys, inners, num_brt_nodes, tree_range, root_level, elements_per_thread](int i) {
		for (int t = i * elements_per_thread+1; t < math::min(num_brt_nodes, (i + 1)*elements_per_thread); ++t)
			MakeNodesHelper(nodes, node_offsets, edge_count, morton_keys, inners, root_level, t, tree_range);
	};


	// Create the threads
	std::vector<std::thread> workers;
	for (int i = 0; i < num_threads; ++i)
		workers.push_back(std::thread(worker_fn, i));
	for (auto& t : workers)	t.join();
}

template<typename DATA_TYPE>
void MakeNodesOpenMP(OctNode<DATA_TYPE>* nodes, const int* node_offsets, const int* edge_count,
               const Code_t* morton_keys, const brt::InnerNodes* inners,
               int num_brt_nodes, int num_threads, float tree_range = 1.0){
    const int root_level = inners[0].delta_node / 3;
    const Code_t root_prefix = morton_keys[0] >> (kCodeLen - (root_level * 3));

    nodes[0].cornor = CodeToPoint(root_prefix << (kCodeLen - (root_level * 3)));
    nodes[0].cell_size = tree_range;

    const auto elements_per_thread = num_brt_nodes / num_threads;
    // skipping root

    #pragma omp parallel num_threads(num_threads)
    {
        int thread_id = omp_get_thread_num();
        int start = thread_id * elements_per_thread;
        int end = (thread_id == num_threads - 1) ? num_brt_nodes : (thread_id + 1) * elements_per_thread;

        #pragma omp for
        for (int t = start; t < end; ++t) {
            MakeNodesHelper(nodes, node_offsets, edge_count, morton_keys, inners, root_level, t, tree_range);
        }
    }
}


template<typename DATA_TYPE>
void MakeNodesHelper(OctNode<DATA_TYPE>* nodes, const int* node_offsets, const int* edge_count,
               const Code_t* morton_keys, const brt::InnerNodes* inners,
               const int root_level, int i, float tree_range = 1.0f){
    int oct_idx = node_offsets[i];
    const int n_new_nodes = edge_count[i];
    for (int j = 0; j < n_new_nodes - 1; ++j) {
      const int level = inners[i].delta_node / 3 - j;
      const Code_t node_prefix = morton_keys[i] >> (kCodeLen - (3 * level));
      const int child_idx = static_cast<int>(node_prefix & 0b111);
      const int parent = oct_idx + 1;

      nodes[parent].SetChild(oct_idx, child_idx);

      // calculate corner point (LSB have already been shifted off)
      nodes[oct_idx].cornor =
          CodeToPoint(node_prefix << (kCodeLen - (3 * level)));

      // each cell is half the size of the level above it
      nodes[oct_idx].cell_size =
          tree_range / static_cast<float>(1 << (level - root_level));

      oct_idx = parent;
    }

    if (n_new_nodes > 0) {
      int rt_parent = inners[i].parent;
      while (edge_count[rt_parent] == 0) {
        rt_parent = inners[rt_parent].parent;
      }
      const int oct_parent = node_offsets[rt_parent];
      const int top_level = inners[i].delta_node / 3 - n_new_nodes + 1;
      const Code_t top_node_prefix =
          morton_keys[i] >> (kCodeLen - (3 * top_level));
      const int child_idx = static_cast<int>(top_node_prefix & 0b111);

      nodes[oct_parent].SetChild(oct_idx, child_idx);
      nodes[oct_idx].cornor =
          CodeToPoint(top_node_prefix << (kCodeLen - (3 * top_level)));
      nodes[oct_idx].cell_size =
          tree_range / static_cast<float>(1 << (top_level - root_level));
    }
}

/**
 * @brief Link the octree nodes together.
 *
 * @param nodes: array of preallocated octree nodes
 * @param node_offsets: ranges of each RT node
 * @param edge_count: number of nodes in each RT node
 * @param morton_keys: sorted morton keys
 * @param inners: binary radix tree nodes
 * @param num_brt_nodes: number of binary radix tree nodes
 */

template<typename DATA_TYPE>
void LinkNodesThreaded(OctNode<DATA_TYPE>* nodes, const int* node_offsets, const int* edge_count,
               const Code_t* morton_keys, const brt::InnerNodes* inners,
               int num_brt_nodes, int num_threads){

	const auto elements_per_thread = math::divide_ceil<int>(num_brt_nodes, num_threads);

  	const auto worker_fn = [nodes, node_offsets, edge_count, morton_keys, inners, num_brt_nodes, elements_per_thread](int i) {
		for (int t = i * elements_per_thread; t < math::min(num_brt_nodes, (i + 1)*elements_per_thread); ++t)
			LinkNodesHelper(nodes, node_offsets, edge_count, morton_keys, inners, t);
	};


	// Create the threads
	std::vector<std::thread> workers;
	for (int i = 0; i < num_threads; ++i)
		workers.push_back(std::thread(worker_fn, i));
	for (auto& t : workers)	t.join();
}

template<typename DATA_TYPE>
void LinkNodesOpenMP(OctNode<DATA_TYPE>* nodes, const int* node_offsets, const int* edge_count,
               const Code_t* morton_keys, const brt::InnerNodes* inners,
               int num_brt_nodes, int num_threads){

    const auto elements_per_thread = num_brt_nodes / num_threads;

    #pragma omp parallel num_threads(num_threads)
    {
        int thread_id = omp_get_thread_num();
        int start = thread_id * elements_per_thread;
        int end = (thread_id == num_threads - 1) ? num_brt_nodes : (thread_id + 1) * elements_per_thread;

        #pragma omp for
        for (int t = start; t < end; ++t) {
            LinkNodesHelper(nodes, node_offsets, edge_count, morton_keys, inners, t);
        }
    }
}

template<typename DATA_TYPE>
void LinkNodesHelper(OctNode<DATA_TYPE>* nodes, const int* node_offsets, const int* edge_count,
               const Code_t* morton_keys, const brt::InnerNodes* inners,
               int i){
    if (IsLeaf(inners[i].left)) {
      const int leaf_idx = GetLeafIndex(inners[i].left);
      const int leaf_level = inners[i].delta_node / 3 + 1;
      const Code_t leaf_prefix =
          morton_keys[leaf_idx] >> (kCodeLen - (3 * leaf_level));

      const int child_idx = static_cast<int>(leaf_prefix & 0b111);
      // walk up the radix tree until finding a node which contributes an
      // octnode
      int rt_node = i;
      while (edge_count[rt_node] == 0) {
        rt_node = inners[rt_node].parent;
      }
      // the lowest octnode in the string contributed by rt_node will be the
      // lowest index
      const int bottom_oct_idx = node_offsets[rt_node];
      nodes[bottom_oct_idx].SetLeaf(leaf_idx, child_idx);
    }

    if (IsLeaf(inners[i].right)) {
      const int leaf_idx = GetLeafIndex(inners[i].left) + 1;
      const int leaf_level = inners[i].delta_node / 3 + 1;
      const Code_t leaf_prefix =
          morton_keys[leaf_idx] >> (kCodeLen - (3 * leaf_level));

      const int child_idx = static_cast<int>(leaf_prefix & 0b111);

      // walk up the radix tree until finding a node which contributes an
      // octnode
      int rt_node = i;
      while (edge_count[rt_node] == 0) {
        rt_node = inners[rt_node].parent;
      }
      // the lowest octnode in the string contributed by rt_node will be the
      // lowest index
      const int bottom_oct_idx = node_offsets[rt_node];
      nodes[bottom_oct_idx].SetLeaf(leaf_idx, child_idx);
    }
}


template<typename DATA_TYPE>
void CheckTree(const Code_t prefix, const int code_len, const OctNode<DATA_TYPE>* nodes,
               const int oct_idx, const Code_t* codes) {
  const OctNode<DATA_TYPE>& node = nodes[oct_idx];
  for (int i = 0; i < 8; ++i) {
    Code_t new_pref = (prefix << 3) | i;
    if (node.child_node_mask & (1 << i)) {
      CheckTree(new_pref, code_len + 3, nodes, node.children[i], codes);
    }
    if (node.child_leaf_mask & (1 << i)) {
      Code_t leaf_prefix =
          codes[node.children[i]] >> (kCodeLen - (code_len + 3));
      if (new_pref != leaf_prefix) {
        printf("oh no...\n");
      }
    }
  }
}

}  // namespace oct