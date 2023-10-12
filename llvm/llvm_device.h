
#pragma once

#include "common/device.h" 
namespace redwood::lang {

class JITModule;
struct LLVMRuntime;

class LlvmDevice : public Device {
 public:
  struct LlvmRuntimeAllocParams : AllocParams {
    JITModule *runtime_jit{nullptr};
    LLVMRuntime *runtime{nullptr};
    uint64_t *result_buffer{nullptr};
    bool use_memory_pool{false};
  };

  template <typename DEVICE>
  DEVICE *as() {
    auto *device = dynamic_cast<DEVICE *>(this);
    RW_ASSERT(device != nullptr);
    return device;
  }

  virtual void *get_memory_addr(DeviceAllocation devalloc) {
    RW_NOT_IMPLEMENTED
  }

  virtual std::size_t get_total_memory() {
    RW_NOT_IMPLEMENTED
  }

  virtual DeviceAllocation import_memory(void *ptr, size_t size) {
    RW_NOT_IMPLEMENTED
  }

  virtual DeviceAllocation allocate_memory_runtime(
      const LlvmRuntimeAllocParams &params) {
    RW_NOT_IMPLEMENTED;
  }

  virtual void clear() {
    RW_NOT_IMPLEMENTED;
  }

  virtual uint64_t *allocate_llvm_runtime_memory_jit(
      const LlvmRuntimeAllocParams &params) {
    RW_NOT_IMPLEMENTED;
  }
};

}  