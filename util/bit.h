#include <cstdint>

namespace redwood{
inline __attribute__((always_inline)) constexpr uint32_t log2int(uint64_t value) {
  int ret = 0;
  value >>= 1;
  while (value) {
    value >>= 1;
    ret += 1;
  }
  return ret;
}

inline __attribute__((always_inline)) constexpr bool is_power_of_two(int32_t x) {
  return x != 0 && (x & (x - 1)) == 0;
}

inline __attribute__((always_inline)) constexpr bool is_power_of_two(uint32_t x) {
  return x != 0 && (x & (x - 1)) == 0;
}

inline __attribute__((always_inline)) constexpr bool is_power_of_two(int64_t x) {
  return x != 0 && (x & (x - 1)) == 0;
}

inline __attribute__((always_inline)) constexpr bool is_power_of_two(uint64_t x) {
  return x != 0 && (x & (x - 1)) == 0;
}
}