#include "llvm/device_memory_pool.h"

#include <memory>


#include "cuda/cuda_driver.h"
#include "cuda/cuda_device.h"
#include <sys/mman.h>

namespace redwood::lang {

DeviceMemoryPool::DeviceMemoryPool(bool merge_upon_release)
    : merge_upon_release_(merge_upon_release) {
  allocator_ = std::make_unique<CachingAllocator>(merge_upon_release);
}

void *DeviceMemoryPool::allocate_with_cache(
    LlvmDevice *device,
    const LlvmDevice::LlvmRuntimeAllocParams &params) {
  std::lock_guard<std::mutex> _(mut_allocation_);

  return allocator_->allocate(device, params);
}

void *DeviceMemoryPool::allocate(std::size_t size,
                                 std::size_t alignment,
                                 bool managed) {
  std::lock_guard<std::mutex> _(mut_allocation_);

  return allocate_raw_memory(size, managed);
}

void DeviceMemoryPool::release(std::size_t size, void *ptr, bool release_raw) {
  std::lock_guard<std::mutex> _(mut_allocation_);

  if (release_raw) {
    deallocate_raw_memory(ptr);
  } else {
    allocator_->release(size, (uint64_t *)ptr);
  }
}

void *DeviceMemoryPool::allocate_raw_memory(std::size_t size, bool managed) {
  /*
    Be aware that this methods is not protected by the mutex.

    allocate_raw_memory() is designed to be a private method, and
    should only be called by its Allocators friends.

    The caller ensures that no other thread is accessing the memory pool
    when calling this method.
  */
  void *ptr = nullptr;

  if (!managed) {
    CUDADriver::get_instance().malloc(&ptr, size);
  } else {
    CUDADriver::get_instance().malloc_managed(&ptr, size, CU_MEM_ATTACH_GLOBAL);
  }


  if (ptr == nullptr) {
    RW_ERROR("Device memory allocation ({} B) failed.", size);
  }

  if (raw_memory_chunks_.count(ptr)) {
    RW_ERROR("Memory address ({:}) is already allocated", ptr);
  }

  raw_memory_chunks_[ptr] = size;
  return ptr;
}

void DeviceMemoryPool::deallocate_raw_memory(void *ptr) {
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

  CUDADriver::get_instance().mem_free(ptr);
  raw_memory_chunks_.erase(ptr);

}

void DeviceMemoryPool::reset() {
  std::lock_guard<std::mutex> _(mut_allocation_);

  const auto ptr_map_copied = raw_memory_chunks_;
  for (auto &ptr : ptr_map_copied) {
    deallocate_raw_memory(ptr.first);
  }
  allocator_ = std::make_unique<CachingAllocator>(merge_upon_release_);
}

DeviceMemoryPool::~DeviceMemoryPool() {
  reset();
}

const size_t DeviceMemoryPool::page_size{1 << 12};  // 4 KB page size by default

DeviceMemoryPool &DeviceMemoryPool::get_instance(bool merge_upon_release) {
  static DeviceMemoryPool *cuda_memory_pool =
      new DeviceMemoryPool(merge_upon_release);
  return *cuda_memory_pool;
}

}  // namespace redwood::lang