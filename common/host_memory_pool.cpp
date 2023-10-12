#include "common/host_memory_pool.h"
#include <memory>

#include <sys/mman.h>


namespace redwood::lang {

HostMemoryPool::HostMemoryPool() {
  allocator_ = std::unique_ptr<UnifiedAllocator>(new UnifiedAllocator());

  RW_TRACE("Memory pool created. Default buffer size per allocator = {} MB",
           UnifiedAllocator::default_allocator_size / 1024 / 1024);
}

void *HostMemoryPool::allocate(std::size_t size,
                               std::size_t alignment,
                               bool exclusive) {
  std::lock_guard<std::mutex> _(mut_allocation_);

  if (!allocator_) {
    RW_ERROR("Memory pool is already destroyed");
  }
  void *ret = allocator_->allocate(size, alignment, exclusive);
  return ret;
}

void HostMemoryPool::release(std::size_t size, void *ptr) {
  std::lock_guard<std::mutex> _(mut_allocation_);

  if (!allocator_) {
    RW_ERROR("Memory pool is already destroyed");
  }

  if (allocator_->release(size, ptr)) {
    if (dynamic_cast<UnifiedAllocator *>(allocator_.get())) {
      deallocate_raw_memory(ptr);  // release raw memory as well
    }
  }
}

void *HostMemoryPool::allocate_raw_memory(std::size_t size) {
  /*
    Be aware that this methods is not protected by the mutex.

    allocate_raw_memory() is designed to be a private method, and
    should only be called by its Allocators friends.

    The caller ensures that no other thread is accessing the memory pool
    when calling this method.
  */

  void *ptr = nullptr;
  ptr = mmap(nullptr, size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS,
             -1, 0);
  RW_ERROR_IF(ptr == MAP_FAILED, "Virtual memory allocation ({} B) failed.",
              size);

  RW_ERROR_IF(((uint64_t)ptr) % page_size != 0,
              "Allocated address ({:}) is not aligned by page size {}", ptr,
              page_size);

  if (raw_memory_chunks_.count(ptr)) {
    RW_ERROR("Memory address ({:}) is already allocated", ptr);
  }

  raw_memory_chunks_[ptr] = size;
  return ptr;
}

void HostMemoryPool::deallocate_raw_memory(void *ptr) {
  /*
    Be aware that this methods is not protected by the mutex.

    deallocate_raw_memory() is designed to be a private method, and
    should only be called by its Allocators friends.

    The caller ensures that no other thread is accessing the memory pool
    when calling this method.
  */
  if (!raw_memory_chunks_.count(ptr)) {
    RW_ERROR("Memory address ({:}) is not allocated", ptr);
  }

  std::size_t size = raw_memory_chunks_[ptr];
  if (munmap(ptr, size) != 0)

    RW_ERROR("Failed to free virtual memory ({} B)", size);

  raw_memory_chunks_.erase(ptr);
}

void HostMemoryPool::reset() {
  std::lock_guard<std::mutex> _(mut_allocation_);
  allocator_ = std::unique_ptr<UnifiedAllocator>(new UnifiedAllocator());

  const auto ptr_map_copied = raw_memory_chunks_;
  for (auto &ptr : ptr_map_copied) {
    deallocate_raw_memory(ptr.first);
  }
}

HostMemoryPool::~HostMemoryPool() {
  reset();
}

const size_t HostMemoryPool::page_size{1 << 12};  // 4 KB page size by default

HostMemoryPool &HostMemoryPool::get_instance() {
  static HostMemoryPool *memory_pool = new HostMemoryPool();
  return *memory_pool;
}

}  // namespace redwood::lang