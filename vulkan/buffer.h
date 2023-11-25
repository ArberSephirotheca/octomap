#pragma once

#include <spdlog/spdlog.h>

#include <numeric>
#include <vulkan/vulkan.h>

#include "vulkan/vma_allocator.h"
#include "vulkan/vulkan_resource.h"
#include "common/logging.h"

namespace redwood::lang {
namespace vulkan {
class Buffer;

using BufferReference = std::reference_wrapper<const Buffer>;

/**
 * @brief Buffer class is an abstraction of data, in the sense of a block of
 * memory that will be processes by the GPU. It is a wrapper of Vulkan buffer
 * managed by Vulkan Memory Allocator (VMA).
 *
 * I provide some default flags for the buffer usage and memory usage. I wanted
 * to use the unified shared memory (USM), for integrated GPUs. The data is
 * shared between CPU and GPU.
 */
class Buffer final : public VulkanResource<vk::Buffer> {
 public:
  /**
   * @brief Construct a new Data Buffer object of any type. By default it is
   * unified shared memory. Make sure to specify the size in total memory size,
   * not data size. (N*sizeof(T) instead of N)
   *
   * @param device_ptr Pointer to the device
   * @param size Size of the buffer. It is the total memory size of the buffer,
   * not the data size.
   * @param buffer_usage Usage of the buffer. Default is storage buffer.
   * @param memory_usage Memory usage of the buffer. Default is Auto.
   * @param flags Allocation flags. Default is persistent mapped memory.
   */
  explicit Buffer(std::shared_ptr<vk::Device> device_ptr,
                  vk::DeviceSize size,
                  vk::BufferUsageFlags buffer_usage =
                      vk::BufferUsageFlagBits::eStorageBuffer,
                  VmaMemoryUsage memory_usage = VMA_MEMORY_USAGE_AUTO,
                  VmaAllocationCreateFlags flags =
                      VMA_ALLOCATION_CREATE_HOST_ACCESS_RANDOM_BIT |
                      VMA_ALLOCATION_CREATE_MAPPED_BIT);

  Buffer(const Buffer &) = delete;

  ~Buffer() override {
    RW_DEBUG("Buffer::~Buffer");
    destroy();
  }

  void destroy() override;

  Buffer &operator=(const Buffer &) = delete;
  Buffer &operator=(Buffer &&) = delete;

  // ---------------------------------------------------------------------------
  //                            Getters
  // ---------------------------------------------------------------------------

  [[nodiscard]] VmaAllocation get_allocation() const { return allocation_; }

  [[nodiscard]] const std::byte *get_data() const { return mapped_data_; }
  [[nodiscard]] std::byte *get_data_mut() { return mapped_data_; }

  template <typename T>
  [[nodiscard]] T *get_data_mut() {
    return reinterpret_cast<T *>(mapped_data_);
  }

  [[nodiscard]] vk::DeviceMemory get_memory() const { return memory_; }
  [[nodiscard, maybe_unused]] vk::DeviceAddress get_device_address() const {
    return device_ptr_->getBufferAddressKHR(get_handle());
  }
  [[nodiscard]] vk::DeviceSize get_size() const { return size_; }

  // ---------------------------------------------------------------------------
  //              Methods for manipulating data in the buffer
  // ---------------------------------------------------------------------------

  template <typename T>
  void tmp_debug_data(const size_t size, const size_t offset = 0) const {
    auto *ptr = reinterpret_cast<T *>(mapped_data_ + offset);
    std::iota(ptr, ptr + size / sizeof(T), T());
  }

  void tmp_write_data(const void *data,
                      const size_t size,
                      const size_t offset = 0) const {
    RW_INFO("Writing {} bytes to buffer", size);
    std::memcpy(mapped_data_ + offset, data, size);
  }

  void tmp_fill_zero(const size_t size, const size_t offset = 0) const {
    RW_INFO("Filling zeros {} bytes to buffer", size);
    std::memset(mapped_data_ + offset, 0, size);
  }

  // ---------------------------------------------------------------------------
  //      The following functions provides infos for the descriptor set
  // ---------------------------------------------------------------------------

  [[nodiscard]] vk::DescriptorBufferInfo construct_descriptor_buffer_info()
      const;

 private:
  // Vulkan Memory Allocator components
  VmaAllocation allocation_ = VK_NULL_HANDLE;
  vk::DeviceMemory memory_ = nullptr;
  vk::DeviceSize size_ = 0;

  // Raw pointer to the mapped data, CPU/GPU shared memory.
  std::byte *mapped_data_ = nullptr;

  bool persistent_ = true;
};
}  // namespace vulkan
}  // namespace redwood::lang