
#include "binary_radix_tree.h"
#include "math.h"
#include <assert.h>     /* assert */
#include <thread>
#include <vector>		/* to store references to thread objects */

namespace binary_radix_tree
{
	//In the context of octree construction or spatial indexing, 
	// this function could be used to calculate the level of the octree at which two points diverge or differ in their Morton code representation.
	/** length of common prefix in bits */
	template <typename T>
	inline int delta(const std::vector<T>& morton_keys, int depth, const int& i, const int& j)
	{
		const int number_of_bits = sizeof(T) * 8;
		//const int number_of_effective_bits = depth*3;
		for (int b = number_of_bits - 1; b >= 0; --b)
			if ((morton_keys[i] >> b) != (morton_keys[j] >> b))
				return static_cast<int>(number_of_bits - b - 1);
		return number_of_bits;
	}

	/** does basic out of bounds checking before invoking delta */
	inline int delta_safe(const int key_num, const std::vector<uint32_t>& morton_keys, int depth, const int& i, const int& j)
	{
		/// for simplicity we define delta(j) to be -1 when j != [0, n)
		return (j < 0 || j >= key_num) ? -1 : delta(morton_keys, depth, i, j);
	}
		
	/** process a single node. is thread safe. */
	void process_internal_node(const int key_num, const std::vector<uint32_t>& morton_keys, int depth, const int i, binary_radix_tree::brt& brt)
	{
		// ________________________________________________________________
		// determine the direction of the range = (+1 or -1)
		int left_delta = delta_safe(key_num,morton_keys, depth, i, i + 1);
		int right_delta = delta_safe(key_num, morton_keys, depth, i, i - 1);

		const auto direction{ math::sign<int>(delta_safe(key_num,morton_keys, depth, i, i + 1) - delta_safe(key_num, morton_keys, depth, i, i - 1)) };
		// sanity check. if it is zero, the set has duplicate data, and we do not resolve duplicate data in this implementation
		assert(direction != 0); 
		
		// ________________________________________________________________
		// compute the upper bound for the length of the edge

		const auto delta_min{ delta_safe(key_num, morton_keys, depth, i, i - direction) };
		int I_max{ 2 };
		while (delta_safe(key_num, morton_keys, depth, i, i + I_max*direction) > delta_min)
			I_max <<= 2;

		// ________________________________________________________________
		// find the other end using binary search

		int I{ 0 };
		for (int t{ I_max / 2 }; t; t /= 2)
			if (delta_safe(key_num, morton_keys, depth, i, i + (I + t)*direction) > delta_min)
				I += t;
		const int j{ i + I * direction };

		// ________________________________________________________________
		// find the split position using binary search

		const auto delta_node{ delta_safe(key_num, morton_keys, depth, i, j) }; // the distance of the prefix of i
		std::cout<<"delta node "<<i <<": "<<delta_node<<std::endl;
		binary_radix_tree::node::set_delta_value(i,brt.nodes_delta, delta_node);
		auto s{ 0 };
		
		int t{ I };
		do
		{
			t = math::divide2ceil<int>(t);
			if (delta_safe(key_num, morton_keys, depth, i, i + (s + t)*direction) > delta_node)
				s += t;

		} while (t > 1);
		
		// split = split position
		const auto split{ i + s * direction + math::min<int>(direction, 0) };

		// ________________________________________________________________
		// sanity check
		assert(delta(morton_keys, depth, i,j) > delta_min);
		assert(delta(morton_keys, depth, split, split + 1) == delta(morton_keys, depth, i, j));
		assert(!(split < 0 || split+1 >= key_num));

		// ________________________________________________________________
		// output child pointers	

		const int left{ (math::min<int>(i, j) == split) ? binary_radix_tree::node::make_leaf<int>(split) : binary_radix_tree::node::make_internal<int>(split) };
		const int right{ (math::max<int>(i, j) == split + 1) ? binary_radix_tree::node::make_leaf<int>(split + 1) : binary_radix_tree::node::make_internal<int>(split + 1) };
		binary_radix_tree::node::set_internal_values<int>(i, brt.internal_nodes, left, right);
	}
};


// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------


/** given a sorted sequence of morton_keys, create a binary radix tree. */
binary_radix_tree::brt binary_radix_tree::create(const int key_num, const std::vector<uint32_t>& morton_keys, int depth)
{
	// initialize the memory to hold the tree
	binary_radix_tree::brt tree(key_num);

	// for each internal node with index i [0, n-2] (to be launched in parallel)
	for (int i = 0; i < key_num-1; ++i)
		process_internal_node(key_num, morton_keys, depth, i, tree);

	return tree;
}

/** given a sorted sequence of morton_keys, create a binary radix tree. */
binary_radix_tree::brt binary_radix_tree::create_threaded(const int key_num, const std::vector<uint32_t>& morton_keys, int depth, const int thread_number)
{
	// initialize the memory to hold the tree
	binary_radix_tree::brt tree(key_num);
	// compute the number of elements a thread will cover
	const auto elements_per_thread{ math::divideceil<int>(key_num, thread_number) };
	
	const auto worker_fn = [key_num, &morton_keys, depth, elements_per_thread, &tree](int i) {
		for (int t = i * elements_per_thread; t < math::min(key_num - 1, (i + 1)*elements_per_thread); ++t)
			process_internal_node(key_num, morton_keys, depth, t, tree);
	};

	// Create the threads
	std::vector<std::thread> workers;
	for (int i = 0; i < thread_number; ++i)
		workers.push_back(std::thread(worker_fn, i));
	for (auto& t : workers)	t.join();	

	return tree;
}

// Function to calculate the count of subprefixes divisible by 3
int countDivisibleSubprefixes(int deltaParent, int deltaChild) {
    int countChild = std::floor(deltaChild / 3);
    int countParent = std::floor(deltaParent / 3);
    return countChild - countParent;
}

// newly added
void binary_radix_tree::create_octree_nodes(const int key_num, const binary_radix_tree::brt& brt){
	std::vector<int> counts(key_num, 0);
	counts[0] = 1;
	int cnt;
	int left_idx;
	int right_idx;
	bool leaf;
	for(int i = 0; i < key_num-1; ++i){
		cnt = 0;
		// get the length of prefix bit of current internal node
		int parnet_prefix = brt.nodes_delta[i];
		std::cout << "prefix of internal node "<< i <<": "<<parnet_prefix<<std::endl;
		// get the left child index at sorted morton code array
		binary_radix_tree::node::get_internal_left(i, brt.internal_nodes, left_idx, leaf);
		std::cout <<"left idx: " <<left_idx<<std::endl;
		int left_prefix = brt.nodes_delta[left_idx];
		std::cout <<"prefix of left child node: "<<left_prefix<<std::endl;
		//  get the right child index at sorted morton code array
		binary_radix_tree::node::get_internal_right(i, brt.internal_nodes, right_idx, leaf);
		int right_prefix = brt.nodes_delta[right_idx];
		std::cout <<"right idx: " <<right_idx<<std::endl;
		std::cout <<"prefix of right child node: "<<right_prefix<<std::endl<<std::endl;
		//cnt += countDivisibleSubprefixes(parnet_prefix, left_prefix);
		//cnt += countDivisibleSubprefixes(parnet_prefix, right_prefix);
		
		// counts ony when (delta_parent-delta_child) is divisible by three
		if((left_prefix-parnet_prefix) %3 == 0){
			for(int j = parnet_prefix+1; j <= left_prefix; ++j){
				if(j%3 == 0){
					cnt += 1;
				}
			}
		}
		if((right_idx-parnet_prefix) %3 == 0){
			for(int j = parnet_prefix+1; j <= right_idx; ++j){
				if(j %3 == 0){
					cnt += 1;
				}
			}
		}		
		counts[i+1] = cnt;
	}

	for(int i = 0; i < key_num; i++){
		std::cout<<"counts[" << i<<"]"<<counts[i]<<std::endl;
	}
}

/** */
unsigned int binary_radix_tree::hardware_concurrency()
{
	return std::thread::hardware_concurrency();
}



// -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------





/** constructor */
binary_radix_tree::brt::brt(const int key_num) :
	internal_nodes(new int[(key_num - 1)*2]),
	leaf_nodes(new int[(key_num-1)*2]), // store the child index of each internal node 
	nodes_delta(new int[(key_num - 1)]),
	key_num(key_num)
{
}


/** move constructor */
binary_radix_tree::brt::brt(brt&& other)
	: internal_nodes(other.internal_nodes)
	, leaf_nodes(other.leaf_nodes)
	, nodes_delta(other.nodes_delta)
	, key_num(other.key_num)
{
	other.internal_nodes = nullptr;
	other.leaf_nodes = nullptr;
	other.nodes_delta = nullptr;
	const_cast<int&>(other.key_num) = 0;
}

/* move assignment operator */
binary_radix_tree::brt & binary_radix_tree::brt::operator=(brt && other)
{
	if (this != &other) {           
		delete[] this->internal_nodes; 
		delete[] this->leaf_nodes;
		delete[] this->nodes_delta;
		this->internal_nodes = other.internal_nodes;
		this->leaf_nodes = other.leaf_nodes;
		this->nodes_delta = other.nodes_delta;
		const_cast<int&>(this->key_num) = other.key_num;
		other.internal_nodes = nullptr;
		other.leaf_nodes = nullptr;
		other.nodes_delta = nullptr;
		const_cast<int&>(other.key_num) = 0;
	}
	return *this;
}
