#pragma once

#include "morton_util.hpp"

namespace brt {

struct InnerNodes {
  // The number of bits in the morton code, this node represents in [Karras]
  int delta_node;

  // pointers
  int left = -1;  // can be either inner or leaf
  int right = -1;
  int parent = -1;
};

// Count leading zeros using bit manipulation
  inline unsigned CountLeadingZeros(const uint64_t num) {
    /*
    if (num == 0) return 64;  // Special case for zero input

    unsigned count = 0;
    uint64_t mask = static_cast<uint64_t>(1) << 63;  // Start from the most significant bit

    while ((num & mask) == 0) {
        count++;
        mask >>= 1;
    }

    return count;
    */
   return __builtin_clzll(num);
}


/**
 * @brief Calculate the number of common prefix bits between two morton codes.
 *
 * @param morton_keys: sorted (not necessary) morton codes
 * @param i: index of the first morton code
 * @param j: index of the second morton code
 * @return number of common prefix bits
 */
 inline int Delta(const Code_t* morton_keys, const int i,
                            const int j) noexcept {
  const auto li = morton_keys[i];
  const auto lj = morton_keys[j];
  constexpr auto unused_bits = 1;
  return CountLeadingZeros(li ^ lj) - unused_bits;
}

/**
 * @brief Calculate the number of common prefix bits between two morton codes.
 * Safe version, return -1 if the index is out of range.
 *
 * @param key_num: Number of morton codes (after remove duplicate)
 * @param morton_keys: sorted (not necessary) morton codes
 * @param i: index of the first morton code
 * @param j: index of the second morton code
 * @return number of common prefix bits
 */
 inline int DeltaSafe(const int key_num, const Code_t* morton_keys,
                                const int i, const int j) noexcept {
  return (j < 0 || j >= key_num) ? -1 : Delta(morton_keys, i, j);
}

/**
 * @brief Given a sorted array of morton codes, make the nodes of binary radix
 * tree. The radix has 'n-1' internal nodes.
 *
 * @param key_num: number of sorted morton codes
 * @param morton_keys: sorted morton codes
 * @param brt_nodes: output an array of internal nodes of size 'n-1'
 */
void create_binary_radix_tree(int key_num, const Code_t* morton_keys,
                          InnerNodes* brt_nodes);

void create_binary_radix_tree_threaded(int key_num, const Code_t* morton_keys,
                          InnerNodes* brt_nodes, int thread_number);
namespace node {

template <typename T>
 T make_leaf(const T& index) {
  return index ^ ((-1 ^ index) & 1UL << (sizeof(T) * 8 - 1));
}

template <typename T>
 T make_internal(const T& index) {
  return index;
}
}  // namespace node

}  // namespace brt

namespace math {
template <typename T>
int sign(T val) {
  return (T(0) < val) - (val < T(0));
}

template <typename T>
T min(const T& x, const T& y) {
  return y ^ ((x ^ y) & -(x < y));
}

template <typename T>
T max(const T& x, const T& y) {
  return x ^ ((x ^ y) & -(x < y));
}

template <typename T>
int divide_ceil(const T& x, const T& y) {
  return (x + y - 1) / y;
}

/** Integer division by two, rounding up */
template <typename T>
int divide2_ceil(const T& x) {
  return (x + 1) >> 1;
}
}  // namespace math