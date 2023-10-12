#pragma once

#include <mutex>
#include <unordered_map>
#include <thread>

//#include "redwood/program/kernel_profiler.h"
#include "cuda/cuda_driver.h"

namespace redwood::lang {

class CUDADriver;

class CUDAContext {
 private:
  void *device_;
  void *context_;
  int dev_count_;
  int compute_capability_;
  std::string mcpu_;
  std::mutex lock_;
  CUDADriver &driver_;
  int max_shared_memory_bytes_;
  bool debug_;
  bool supports_mem_pool_;

 public:
  CUDAContext();

  std::size_t get_total_memory();
  std::size_t get_free_memory();
  std::string get_device_name();

  bool detected() const {
    return dev_count_ != 0;
  }

  void launch(void *func,
              const std::string &task_name,
              std::vector<void *> arg_pointers,
              std::vector<int> arg_sizes,
              unsigned grid_dim,
              unsigned block_dim,
              std::size_t dynamic_shared_mem_bytes);

/*
  void set_profiler(KernelProfilerBase *profiler) {
    profiler_ = profiler;
  }
  */

  void set_debug(bool debug) {
    debug_ = debug;
  }

  std::string get_mcpu() const {
    return mcpu_;
  }

  void *get_context() {
    return context_;
  }

  void make_current() {
    driver_.context_set_current(context_);
  }

  int get_compute_capability() const {
    return compute_capability_;
  }

  bool supports_mem_pool() const {
    return supports_mem_pool_;
  }

  ~CUDAContext();

  class ContextGuard {
   private:
    void *old_ctx_;
    void *new_ctx_;

   public:
    explicit ContextGuard(CUDAContext *new_ctx)
        : old_ctx_(nullptr), new_ctx_(new_ctx->context_) {
      CUDADriver::get_instance().context_get_current(&old_ctx_);
      if (old_ctx_ != new_ctx_)
        new_ctx->make_current();
    }

    ~ContextGuard() {
      if (old_ctx_ != new_ctx_) {
        CUDADriver::get_instance().context_set_current(old_ctx_);
      }
    }
  };

  ContextGuard get_guard() {
    return ContextGuard(this);
  }

  std::unique_lock<std::mutex> get_lock_guard() {
    return std::unique_lock<std::mutex>(lock_);
  }

  static CUDAContext &get_instance();
};

}  // namespace redwood::lang