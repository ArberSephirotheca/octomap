#pragma once

#include <spdlog/spdlog.h>

#include <vulkan/vulkan.hpp>

#include <VkBootstrap.h>

namespace redwood::lang {
    namespace vulkan{

// Need a global allocator for VMA
// extern VmaAllocator g_allocator;

/**
 * @brief Basically do the initializations, save you a lot of time. BaseEngine
 * will setup the Vulkan instance, physical device, logical device etc. For
 * compute shader usage only. By default it will pick an integrated GPU, with
 * compute queue. It should only pick integrated GPU.
 */
class BaseEngine {
 public:
  BaseEngine();

  ~BaseEngine() {
    spdlog::debug("BaseEngine::~BaseEngine");
    destroy();
  }

  void destroy() const;

  // ---------------------------------------------------------------------------
  //                  Getter and Setter
  // ---------------------------------------------------------------------------

  [[nodiscard, maybe_unused]] const vkb::Instance &get_instance() const {
    return instance_;
  }

  [[nodiscard, maybe_unused]] const vkb::Device &get_device() const {
    return device_;
  }

  /**
   * @brief Make a shared pointer to the vk::Device, NOT the vkb::Device.
   *
   * @return std::shared_ptr<vk::Device> A shared pointer to the vk::Device.
   */
  [[nodiscard]] auto get_device_ptr() {
    return std::make_shared<vk::Device>(device_.device);
  }

  [[nodiscard, maybe_unused]] vk::Queue &get_queue() { return queue_; }
  [[nodiscard, maybe_unused]] const vk::Queue &get_queue() const {
    return queue_;
  }

 private:
  void device_initialization();
  void get_queues();
  void vma_initialization() const;

 protected:
  vkb::Instance instance_;
  vkb::Device device_;
  vk::Queue queue_;
};
    }
}  // namespace core