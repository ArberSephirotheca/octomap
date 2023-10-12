#pragma once
#include "common/unified_allocator.h"
#include "common/logging.h"
#include <mutex>
#include <vector>
#include <memory>
#include <thread>

namespace redwood::lang {

// A memory pool that runs on the host

class HostMemoryPool {
 public:
  static const size_t page_size;

  static HostMemoryPool &get_instance();

  void *allocate(std::size_t size,
                 std::size_t alignment,
                 bool exclusive = false);
  void release(std::size_t size, void *ptr);
  void reset();
  HostMemoryPool();
  ~HostMemoryPool();

 protected:
  void *allocate_raw_memory(std::size_t size);
  void deallocate_raw_memory(void *ptr);

  // All the raw memory allocated from OS/Driver
  // We need to keep track of them to guarantee that they are freed
  std::map<void *, std::size_t> raw_memory_chunks_;

  std::unique_ptr<UnifiedAllocator> allocator_;
  std::mutex mut_allocation_;

  friend class UnifiedAllocator;
};

}  // namespace redwood::lang