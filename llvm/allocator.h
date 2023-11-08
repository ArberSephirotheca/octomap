#pragma once

#include "common/core.h"
#include "common/arithmetic.h"
#include "llvm/llvm_device.h"
#include "inc/constants.h"
#include <stdint.h>
#include <map>
#include <set>

namespace redwood::lang {

class CachingAllocator {
 public:
  explicit CachingAllocator(bool merge_upon_release = true);

  uint64_t *allocate(LlvmDevice *device,
                     const LlvmDevice::LlvmRuntimeAllocParams &params);
  void release(size_t sz, uint64_t *ptr);

 private:
  void merge_and_insert(uint8_t *ptr, std::size_t size);

  std::set<std::pair<std::size_t, uint8_t *>> mem_blocks_;
  std::map<uint8_t *, std::size_t> ptr_map_;

  // Allocator options
  bool merge_upon_release_ = true;
};

}  // namespace redwood::lang