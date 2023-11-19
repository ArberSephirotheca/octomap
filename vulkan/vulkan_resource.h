#pragma once

#include <memory>
#include <vulkan/vulkan.hpp>

namespace redwood::lang {
  namespace vulkan{

/**
 * @brief This abstract class provides a pointer to the vk::device, so each
 * Vulkan Components that inherits from this class can access the device. Need
 * to pass in a shared_ptr<vk::Device> to the constructor.
 *
 * @tparam HandleT The Vulkan handle type (e.g., vk::Buffer, vk::Image, etc.)
 */
template <typename HandleT>
class VulkanResource {
 public:
  VulkanResource() = delete;

  explicit VulkanResource(std::shared_ptr<vk::Device> device_ptr)
      : device_ptr_{std::move(device_ptr)} {}

  VulkanResource(const VulkanResource &) = delete;
  VulkanResource(VulkanResource &&) = delete;

  virtual ~VulkanResource() = default;

  VulkanResource &operator=(const VulkanResource &) = delete;
  VulkanResource &operator=(VulkanResource &&) = delete;

  [[nodiscard]] HandleT &get_handle() { return handle_; }
  [[nodiscard]] const HandleT &get_handle() const { return handle_; }

 protected:
  virtual void destroy() = 0;

  std::shared_ptr<vk::Device> device_ptr_;
  HandleT handle_;
};
  }

}  // namespace redwood::lang