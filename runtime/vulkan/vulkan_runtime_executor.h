#pragma once

#include <cstddef>
#include <memory>
#include <stack>

#include "vulkan/vulkan_device.h"
#include "runtime/vulkan/vulkan_offline_cache.h"
#include "runtime/vulkan/snode_tree_buffer_manager.h"
//#include "runtime/vulkan/vulkan_context.h"
#include "struct/snode_tree.h"
#include "program/compile_config.h"
#include "system/threading.h"
#include "runtime/vulkan/runtime.h"
#include "runtime/vulkan/vulkan_struct_compiler.h"
#include "cuda/cuda_device.h"
#include "common/host_memory_pool.h"
#include "cpu/cpu_device.h"
#include "cuda/cuda_driver.h"
#include "vulkan/device_memory_pool.h"
#include "cuda/cuda_context.h"

#define RW_RUNRWME_HOST
#include "program/context.h"
#undef RW_RUNRWME_HOST

namespace redwood::lang {

namespace cuda {
class CudaDevice;
}  // namespace cuda

namespace amdgpu {
class AmdgpuDevice;
}  // namespace amdgpu

namespace cpu {
class CpuDevice;
}  // namespace cpu



class vulkanRuntimeExecutor {
 public:
  vulkanRuntimeExecutor(CompileConfig &config/*, KernelProfilerBase *profiler*/);
  virtual ~vulkanRuntimeExecutor();
  /**
   * Initializes the runtime system for vulkan based backends.
   */
  void materialize_runtime(/*KernelProfilerBase *profiler,*/
                           uint64_t **result_buffer_ptr);

  // SNodeTree Allocation
  void initialize_vulkan_runtime_snodes(
      const vulkanOfflineCache::FieldCacheData &field_cache_data,
      uint64_t *result_buffer);
  SNodeTree* add_snode_tree(std::unique_ptr<SNode> root);
  // Ndarray and ArgPack Allocation
  DeviceAllocation allocate_memory_on_device(std::size_t alloc_size,
                                             uint64_t *result_buffer);

  void deallocate_memory_on_device(DeviceAllocation handle);

  //void check_runtime_error(uint64 *result_buffer);

  uint64_t *get_device_alloc_info_ptr(const DeviceAllocation &alloc);

  const CompileConfig &get_config() const {
    return config_;
  }

  //redwoodvulkanContext *get_vulkan_context();

  vulkanRuntime *get_vulkan_runtime();

  Device *get_compute_device();

  vulkanDevice *vulkan_device();

  void synchronize();

  bool use_device_memory_pool() {
    return use_device_memory_pool_;
  }

 private:
  /* ----------------------- */
  /* ------ Allocation ----- */
  /* ----------------------- */
  int allocate_snode_tree_id();

  template <typename T>
  T fetch_result(int i, uint64_t *result_buffer) {
    return redwood_union_cast_with_different_sizes<T>(
        fetch_result_uint64_t(i, result_buffer));
  }

  template <typename T>
  T fetch_result(char *result_buffer, int offset) {
    return *(T *)(result_buffer + offset);
  }

  DevicePtr get_snode_tree_device_ptr(int tree_id);

  void fill_ndarray(const DeviceAllocation &alloc,
                    std::size_t size,
                    uint32_t data);

  void *preallocate_memory(std::size_t prealloc_size,
                           DeviceAllocationUnique &devalloc);
  void preallocate_runtime_memory();

  /* ------------------------- */
  /* ---- Runtime Helpers ---- */
  /* ------------------------- */
  
  //void print_list_manager_info(void *list_manager, uint64 *result_buffer);
  /*
  void print_memory_profiler_info(
      std::vector<std::unique_ptr<SNodeTree>> &snode_trees_,
      uint64 *result_buffer);
    */
   /*
  template <typename T, typename... Args>
  T runtime_query(const std::string &key,
                  uint64 *result_buffer,
                  Args &&...args) {
    RW_ASSERT(arch_uses_vulkan(config_.arch));

    auto   = get_runtime_jit_module();
    runtime->call<void *>("runtime_" + key, vulkan_runtime_,
                          std::forward<Args>(args)...);
    return redwood_union_cast_with_different_sizes<T>(fetch_result_uint64(
        redwood_result_buffer_runtime_query_id, result_buffer));
  }
  */
  void cache_field(int snode_tree_id, int root_id, const StructCompiler& struct_compiler);
  /* -------------------------- */
  /* ------ Member Access ----- */
  /* -------------------------- */
  void finalize();

  uint64_t fetch_result_uint64_t(int i, uint64_t *result_buffer);
  void destroy_snode_tree(SNodeTree *snode_tree);
  /*
  std::size_t get_snode_num_dynamically_allocated(SNode *snode,
                                                  uint64 *result_buffer);
  */
  //void init_runtime_jit_module(std::unique_ptr<vulkan::Module> module);

 private:
  CompileConfig &config_;

  std::vector<std::unique_ptr<SNodeTree>> snode_trees_;
  std::stack<int> free_snode_tree_ids_;
  std::unique_ptr<vulkanOfflineCache> cache_data_;
  //temporary use
  uint64_t *result_buffer_{nullptr};
  //std::unique_ptr<redwoodvulkanContext> vulkan_context_{nullptr};
  //std::unique_ptr<JITSession> jit_session_{nullptr};
  //JITModule *runtime_jit_module_{nullptr};
  void *vulkan_runtime_{nullptr};

  std::unique_ptr<ThreadPool> thread_pool_{nullptr};
  std::shared_ptr<Device> device_{nullptr};

  std::unique_ptr<SNodeTreeBufferManager> snode_tree_buffer_manager_{nullptr};
  std::unordered_map<int, DeviceAllocation> snode_tree_allocs_;
  DeviceAllocationUnique preallocated_runtime_objects_allocs_ = nullptr;
  DeviceAllocationUnique preallocated_runtime_memory_allocs_ = nullptr;
  std::unordered_map<DeviceAllocationId, DeviceAllocation>
      allocated_runtime_memory_allocs_;

  // good buddy
  // friend vulkanProgramImpl;
  friend SNodeTreeBufferManager;

  bool use_device_memory_pool_ = false;
  bool finalized_{false};
  void materialize_snode_tree(SNodeTree *tree,
                                      uint64_t *result_buffer_ptr);
  //KernelProfilerBase *profiler_ = nullptr;
};

}  // namespace redwood::lang

