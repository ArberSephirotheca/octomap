#include "runtime/llvm/llvm_runtime_executor.h"

#include "common/host_memory_pool.h"
#include "cpu/cpu_device.h"
#include "cuda/cuda_device.h"
#include "cuda/cuda_driver.h"
#include "llvm/device_memory_pool.h"
#include "cuda/cuda_context.h"

namespace redwood::lang {
namespace {
void assert_failed_host(const char *msg) {
  RW_ERROR("Assertion failure: {}", msg);
}

void *host_allocate_aligned(HostMemoryPool *memory_pool,
                            std::size_t size,
                            std::size_t alignment) {
  return memory_pool->allocate(size, alignment);
}

}  // namespace


LlvmRuntimeExecutor::LlvmRuntimeExecutor(CompileConfig &config
                                         /*, KernelProfilerBase **/)
    : config_(config) {
  if (config.arch == Arch::cuda) {
    /*
    if (!is_cuda_api_available()) {
      RW_WARN("No CUDA driver API detected.");
      config.arch = host_arch();
    
    }*/ if (!CUDAContext::get_instance().detected()) {
      RW_WARN("No CUDA device detected.");
      config.arch = host_arch();
    } else {
      // CUDA runtime created successfully
      use_device_memory_pool_ = CUDAContext::get_instance().supports_mem_pool();
    }
  }
  /*
  if (config.kernel_profiler) {
    profiler_ = profiler;
  }
  */

  snode_tree_buffer_manager_ = std::make_unique<SNodeTreeBufferManager>(this);
  thread_pool_ = std::make_unique<ThreadPool>(config.cpu_max_num_threads);

  llvm_runtime_ = nullptr;
  cache_data_ = std::make_unique<LlvmOfflineCache>();
  if (arch_is_cpu(config.arch)) {
    config.max_block_dim = 1024;
    device_ = std::make_shared<cpu::CpuDevice>();

  }
  else if (config.arch == Arch::cuda) {
    int num_SMs{1};
    CUDADriver::get_instance().device_get_attribute(
        &num_SMs, CU_DEVICE_ATTRIBUTE_MULTIPROCESSOR_COUNT, nullptr); 
    int query_max_block_dim{1024};
    CUDADriver::get_instance().device_get_attribute(
        &query_max_block_dim, CU_DEVICE_ATTRIBUTE_MAX_BLOCK_DIM_X, nullptr);
    int version{0};
    CUDADriver::get_instance().driver_get_version(&version);
    int query_max_block_per_sm{16};
    if (version >= 11000) {
      // query this attribute only when CUDA version is above 11.0
      CUDADriver::get_instance().device_get_attribute(
          &query_max_block_per_sm,
          CU_DEVICE_ATTRIBUTE_MAX_BLOCKS_PER_MULTIPROCESSOR, nullptr);
    }

    if (config.max_block_dim == 0) {
      config.max_block_dim = query_max_block_dim;
    }

    if (config.saturating_grid_dim == 0) {
      if (version >= 11000) {
        RW_TRACE("CUDA max blocks per SM = {}", query_max_block_per_sm);
      }
      config.saturating_grid_dim = num_SMs * query_max_block_per_sm * 2;
    }
    /*
    if (config.kernel_profiler) {
      CUDAContext::get_instance().set_profiler(profiler);
    } else {
      CUDAContext::get_instance().set_profiler(nullptr);
    }
    
    CUDAContext::get_instance().set_debug(config.debug);
    */

    if (config.cuda_stack_limit != 0) {
      CUDADriver::get_instance().context_set_limit(CU_LIMIT_STACK_SIZE,
                                                   config.cuda_stack_limit);
    }
    device_ = std::make_shared<cuda::CudaDevice>();
  }
  else {
    RW_NOT_IMPLEMENTED
  }
  /*
  llvm_context_ = std::make_unique<redwoodLLVMContext>(
      config_, arch_is_cpu(config.arch) ? host_arch() : config.arch);
      */
  //jit_session_ = JITSession::create(llvm_context_.get(), config, config.arch);
  //init_runtime_jit_module(llvm_context_->clone_runtime_module());
}

/*
redwoodLLVMContext *LlvmRuntimeExecutor::get_llvm_context() {
  return llvm_context_.get();
}
*/
/*
JITModule *LlvmRuntimeExecutor::create_jit_module(
    std::unique_ptr<llvm::Module> module) {
  return jit_session_->add_module(std::move(module));
}
*/

/*
JITModule *LlvmRuntimeExecutor::get_runtime_jit_module() {
  return runtime_jit_module_;
}
*/

/*
void LlvmRuntimeExecutor::print_list_manager_info(void *list_manager,
                                                  uint64_t *result_buffer) {
  auto list_manager_len = runtime_query<int32>("ListManager_get_num_elements",
                                               result_buffer, list_manager);

  auto element_size = runtime_query<int32>("ListManager_get_element_size",
                                           result_buffer, list_manager);

  auto elements_per_chunk =
      runtime_query<int32>("ListManager_get_max_num_elements_per_chunk",
                           result_buffer, list_manager);

  auto num_active_chunks = runtime_query<int32>(
      "ListManager_get_num_active_chunks", result_buffer, list_manager);

  auto size_MB = 1e-6f * num_active_chunks * elements_per_chunk * element_size;

  fmt::print(
      " length={:n}     {:n} chunks x [{:n} x {:n} B]  total={:.4f} MB\n",
      list_manager_len, num_active_chunks, elements_per_chunk, element_size,
      size_MB);
}
*/

void LlvmRuntimeExecutor::synchronize() {
  if (config_.arch == Arch::cuda) {
    CUDADriver::get_instance().stream_synchronize(nullptr);
  }
  fflush(stdout);
}

uint64_t LlvmRuntimeExecutor::fetch_result_uint64_t(int i, uint64_t *result_buffer) {
  // TODO: We are likely doing more synchronization than necessary. Simplify the
  // sync logic when we fetch the result.
  synchronize();
  uint64_t ret;
  if (config_.arch == Arch::cuda) {
    CUDADriver::get_instance().memcpy_device_to_host(&ret, result_buffer + i,
                                                     sizeof(uint64_t));
  }  else {
    ret = result_buffer[i];
  }
  return ret;
}

/*
std::size_t LlvmRuntimeExecutor::get_snode_num_dynamically_allocated(
    SNode *snode,
    uint64_t *result_buffer) {
  RW_ASSERT(arch_uses_llvm(config_.arch));

  auto node_allocator = 
      runtime_query<void *>("LLVMRuntime_get_node_allocators", result_buffer,
                            llvm_runtime_, snode->id);
  auto data_list = runtime_query<void *>("NodeManager_get_data_list",
                                         result_buffer, node_allocator);

  return (std::size_t)runtime_query<int32>("ListManager_get_num_elements",
                                           result_buffer, data_list);
}
*/

/*
void LlvmRuntimeExecutor::check_runtime_error(uint64_t *result_buffer) {
  synchronize();
  auto *runtime_jit_module = get_runtime_jit_module();
  runtime_jit_module->call<void *>("runtime_retrieve_and_reset_error_code",
                                   llvm_runtime_);
  auto error_code =
      fetch_result<int64>(redwood_result_buffer_error_id, result_buffer);

  if (error_code) {
    std::string error_message_template;

    // Here we fetch the error_message_template char by char.
    // This is not efficient, but fortunately we only need to do this when an
    // assertion fails. Note that we may not have unified memory here, so using
    // "fetch_result" that works across device/host memory is necessary.
    for (int i = 0;; i++) {
      runtime_jit_module->call<void *>("runtime_retrieve_error_message",
                                       llvm_runtime_, i);
      auto c = fetch_result<char>(redwood_result_buffer_error_id, result_buffer);
      error_message_template += c;
      if (c == '\0') {
        break;
      }
    }

    if (error_code == 1) {
      const auto error_message_formatted = format_error_message(
          error_message_template,
          [runtime_jit_module, result_buffer, this](int argument_id) {
            runtime_jit_module->call<void *>(
                "runtime_retrieve_error_message_argument", llvm_runtime_,
                argument_id);
            return fetch_result<uint64_t>(redwood_result_buffer_error_id,
                                        result_buffer);
          });
      throw redwoodAssertionError(error_message_formatted);
    } else {
      RW_NOT_IMPLEMENTED
    }
  }
}
*/

/*
void LlvmRuntimeExecutor::print_memory_profiler_info(
    std::vector<std::unique_ptr<SNodeTree>> &snode_trees_,
    uint64_t *result_buffer) {
  RW_ASSERT(arch_uses_llvm(config_.arch));

  fmt::print("\n[Memory Profiler]\n");

  std::locale::global(std::locale("en_US.UTF-8"));
  // So that thousand separators are added to "{:n}" slots in fmtlib.
  // E.g., 10000 is printed as "10,000".
  // TODO: is there a way to set locale only locally in this function?

  std::function<void(SNode *, int)> visit = [&](SNode *snode, int depth) {
    auto element_list =
        runtime_query<void *>("LLVMRuntime_get_element_lists", result_buffer,
                              llvm_runtime_, snode->id);

    if (snode->type != SNodeType::place) {
      fmt::print("SNode {:10}\n", snode->get_node_type_name_hinted());

      if (element_list) {
        fmt::print("  active element list:");
        print_list_manager_info(element_list, result_buffer);

        auto node_allocator =
            runtime_query<void *>("LLVMRuntime_get_node_allocators",
                                  result_buffer, llvm_runtime_, snode->id);

        if (node_allocator) {
          auto free_list = runtime_query<void *>("NodeManager_get_free_list",
                                                 result_buffer, node_allocator);
          auto recycled_list = runtime_query<void *>(
              "NodeManager_get_recycled_list", result_buffer, node_allocator);

          auto free_list_len = runtime_query<int32>(
              "ListManager_get_num_elements", result_buffer, free_list);

          auto recycled_list_len = runtime_query<int32>(
              "ListManager_get_num_elements", result_buffer, recycled_list);

          auto free_list_used = runtime_query<int32>(
              "NodeManager_get_free_list_used", result_buffer, node_allocator);

          auto data_list = runtime_query<void *>("NodeManager_get_data_list",
                                                 result_buffer, node_allocator);
          fmt::print("  data list:          ");
          print_list_manager_info(data_list, result_buffer);

          fmt::print(
              "  Allocated elements={:n}; free list length={:n}; recycled list "
              "length={:n}\n",
              free_list_used, free_list_len, recycled_list_len);
        }
      }
    }
    for (const auto &ch : snode->ch) {
      visit(ch.get(), depth + 1);
    }
  };

  for (auto &a : snode_trees_) {
    visit(a->root(), 0);
  }
  

  auto total_requested_memory = runtime_query<std::size_t>(
      "LLVMRuntime_get_total_requested_memory", result_buffer, llvm_runtime_);

  fmt::print(
      "Total requested dynamic memory (excluding alignment padding): {:n} B\n",
      total_requested_memory);
}
*/

DevicePtr LlvmRuntimeExecutor::get_snode_tree_device_ptr(int tree_id) {
  DeviceAllocation tree_alloc = snode_tree_allocs_[tree_id];
  return tree_alloc.get_ptr();
}

void LlvmRuntimeExecutor::initialize_llvm_runtime_snodes(
    const LlvmOfflineCache::FieldCacheData &field_cache_data,
    uint64_t *result_buffer) {
  //auto *const runtime_jit = get_runtime_jit_module();
  // By the time this creator is called, "this" is already destroyed.
  // Therefore it is necessary to capture members by values.
  size_t root_size = field_cache_data.root_size;
  const auto snode_metas = field_cache_data.snode_metas;
  const int tree_id = field_cache_data.tree_id;
  const int root_id = field_cache_data.root_id;

  bool all_dense = config_.demote_dense_struct_fors;
  
  for (size_t i = 0; i < snode_metas.size(); i++) {
    if (snode_metas[i].type != SNodeType::dense &&
        snode_metas[i].type != SNodeType::place &&
        snode_metas[i].type != SNodeType::root) {
      all_dense = false;
      break;
    }
    
  }
  

  if (config_.arch == Arch::cuda && use_device_memory_pool() && !all_dense) {
    preallocate_runtime_memory();
  }

  RW_TRACE("Allocating data structure of size {} bytes", root_size);
  std::size_t rounded_size = redwood::iroundup(root_size, redwood_page_size);

  Ptr root_buffer = snode_tree_buffer_manager_->allocate(rounded_size, tree_id,
                                                         result_buffer);
  if (config_.arch == Arch::cuda) {
    CUDADriver::get_instance().memset(root_buffer, 0, rounded_size);

  }  else {
    std::memset(root_buffer, 0, rounded_size);
  }
  auto llvm_runtime = static_cast<redwood::runtime::LLVMRuntime*>(llvm_runtime_);
  DeviceAllocation alloc =
      llvm_device()->import_memory(root_buffer, rounded_size);
  snode_tree_allocs_[tree_id] = alloc;
  redwood::runtime::runtime_initialize_snodes(llvm_runtime, root_size, root_id, (int) snode_metas.size(),
                            tree_id, rounded_size, root_buffer, all_dense);

  for (size_t i = 0; i < snode_metas.size(); i++) {
    if (is_gc_able(snode_metas[i].type)) {
      const auto snode_id = snode_metas[i].id;
      std::size_t node_size;
      auto element_size = snode_metas[i].cell_size_bytes;
      
      if (snode_metas[i].type == SNodeType::pointer) {
        // pointer. Allocators are for single elements
        node_size = element_size;
      }
      // else {
        // dynamic. Allocators are for the chunks
        node_size = sizeof(void *) + element_size * snode_metas[i].chunk_size;
      //}
      RW_TRACE("Initializing allocator for snode {} (node size {})", snode_id,
               node_size);
      redwood::runtime::runtime_NodeAllocator_initialize(llvm_runtime, snode_id, node_size);
      RW_TRACE("Allocating ambient element for snode {} (node size {})",
               snode_id, node_size);
      redwood::runtime::runtime_allocate_ambient(llvm_runtime, snode_id, node_size);
    }
  }
}

LlvmDevice *LlvmRuntimeExecutor::llvm_device() {
  RW_ASSERT(dynamic_cast<LlvmDevice *>(device_.get()));
  return static_cast<LlvmDevice *>(device_.get());
}

DeviceAllocation LlvmRuntimeExecutor::allocate_memory_on_device(
    std::size_t alloc_size,
    uint64_t *result_buffer) {
    LlvmDevice::LlvmRuntimeAllocParams param = LlvmDevice::LlvmRuntimeAllocParams{{alloc_size, /*host_write=*/false, /*host_read=*/false,
        /*export_sharing=*/false, AllocUsage::Storage},
       //get_runtime_jit_module(),
       get_llvm_runtime(),
       result_buffer,
       use_device_memory_pool()};
  auto devalloc = llvm_device()->allocate_memory_runtime(param);

  RW_ASSERT(allocated_runtime_memory_allocs_.find(devalloc.alloc_id) ==
            allocated_runtime_memory_allocs_.end());
  allocated_runtime_memory_allocs_[devalloc.alloc_id] = devalloc;
  return devalloc;
}

void LlvmRuntimeExecutor::deallocate_memory_on_device(DeviceAllocation handle) {
  RW_ASSERT(allocated_runtime_memory_allocs_.find(handle.alloc_id) !=
            allocated_runtime_memory_allocs_.end());
  llvm_device()->dealloc_memory(handle);
  allocated_runtime_memory_allocs_.erase(handle.alloc_id);
}

void LlvmRuntimeExecutor::fill_ndarray(const DeviceAllocation &alloc,
                                       std::size_t size,
                                       uint32_t data) {
  auto ptr = get_device_alloc_info_ptr(alloc);
  if (config_.arch == Arch::cuda) {
    CUDADriver::get_instance().memsetd32((void *)ptr, data, size);

  } else {
    std::fill((uint32_t *)ptr, (uint32_t *)ptr + size, data);
  }
}

uint64_t *LlvmRuntimeExecutor::get_device_alloc_info_ptr(
    const DeviceAllocation &alloc) {
  if (config_.arch == Arch::cuda) {
    return (uint64_t *)llvm_device()
        ->as<cuda::CudaDevice>()
        ->get_alloc_info(alloc)
        .ptr;
  }

  return (uint64_t *)llvm_device()
      ->as<cpu::CpuDevice>()
      ->get_alloc_info(alloc)
      .ptr;
}

void LlvmRuntimeExecutor::finalize() {
  //profiler_ = nullptr;
  if (config_.arch == Arch::cuda || config_.arch == Arch::amdgpu) {
    preallocated_runtime_objects_allocs_.reset();
    preallocated_runtime_memory_allocs_.reset();

    // Reset runtime memory
    auto allocated_runtime_memory_allocs_copy =
        allocated_runtime_memory_allocs_;
    for (auto &iter : allocated_runtime_memory_allocs_copy) {
      // The runtime allocation may have already been freed upon explicit
      // Ndarray/Field destruction Check if the allocation still alive
      void *ptr = llvm_device()->get_memory_addr(iter.second);
      if (ptr == nullptr)
        continue;

      deallocate_memory_on_device(iter.second);
    }
    allocated_runtime_memory_allocs_.clear();

    // Reset device
    llvm_device()->clear();

    // Reset memory pool
    DeviceMemoryPool::get_instance().reset();

    // Release unused memory from cuda memory pool
    synchronize();
  }
  finalized_ = true;
}

LlvmRuntimeExecutor::~LlvmRuntimeExecutor() {
  if (!finalized_) {
    finalize();
  }
}

void *LlvmRuntimeExecutor::preallocate_memory(
    std::size_t prealloc_size,
    DeviceAllocationUnique &devalloc) {
  DeviceAllocation preallocated_device_buffer_alloc;

  Device::AllocParams preallocated_device_buffer_alloc_params;
  preallocated_device_buffer_alloc_params.size = prealloc_size;
  RedwoodResult res =
      llvm_device()->allocate_memory(preallocated_device_buffer_alloc_params,
                                     &preallocated_device_buffer_alloc);
  RW_ERROR_IF(res != RedwoodResult::success,
              "Failed to pre-allocate device memory (err: {})", int(res));

  void *preallocated_device_buffer =
      llvm_device()->get_memory_addr(preallocated_device_buffer_alloc);
  devalloc = std::make_unique<DeviceAllocationGuard>(
      std::move(preallocated_device_buffer_alloc));
  return preallocated_device_buffer;
}


void LlvmRuntimeExecutor::preallocate_runtime_memory() {
  if (preallocated_runtime_memory_allocs_ != nullptr)
    return;

  std::size_t total_prealloc_size = 0;
  const auto total_mem = llvm_device()->get_total_memory();
  if (config_.device_memory_fraction == 0) {
    RW_ASSERT(config_.device_memory_GB > 0);
    total_prealloc_size = std::size_t(config_.device_memory_GB * (1UL << 30));
  } else {
    total_prealloc_size =
        std::size_t(config_.device_memory_fraction * total_mem);
  }
  RW_ASSERT(total_prealloc_size <= total_mem);

  void *runtime_memory_prealloc_buffer = preallocate_memory(
      total_prealloc_size, preallocated_runtime_memory_allocs_);

  RW_TRACE("Allocating device memory {:.2f} MB",
           1.0 * total_prealloc_size / (1UL << 20));

  auto llvm_runtime = static_cast<redwood::runtime::LLVMRuntime*>(llvm_runtime_);

  redwood::runtime::runtime_initialize_memory(llvm_runtime, total_prealloc_size,
                              (Ptr) runtime_memory_prealloc_buffer);
  /*
  auto *const runtime_jit = get_runtime_jit_module();
  runtime_jit->call<void *, std::size_t, void *>(
      "runtime_initialize_memory", llvm_runtime_, total_prealloc_size,
      runtime_memory_prealloc_buffer);
      */
}


void LlvmRuntimeExecutor::materialize_runtime(/*KernelProfilerBase *profiler,*/
                                              uint64_t **result_buffer_ptr) {
  // Starting random state for the program calculated using the random seed.
  // The seed is multiplied by 1048391 so that two programs with different seeds
  // will not have overlapping random states in any thread.
  int starting_rand_state = config_.random_seed * 1048391;

  // Number of random states. One per CPU/CUDA thread.
  int num_rand_states = 0;

  if (config_.arch == Arch::cuda || config_.arch == Arch::amdgpu) {
#if defined(RW_WITH_CUDA) || defined(RW_WITH_AMDGPU)
    // It is important to make sure that every CUDA thread has its own random
    // state so that we do not need expensive per-state locks.
    num_rand_states = config_.saturating_grid_dim * config_.max_block_dim;
#else
    RW_NOT_IMPLEMENTED
#endif
  } else {
    num_rand_states = config_.cpu_max_num_threads;
  }

  // The result buffer allocated here is only used for the launches of
  // runtime JIT functions. To avoid memory leak, we use the head of
  // the preallocated device buffer as the result buffer in
  // CUDA and AMDGPU backends.
  // | ==================preallocated device buffer ========================== |
  // |<- reserved for return ->|<---- usable for allocators on the device ---->|
  //auto *const runtime_jit = get_runtime_jit_module();

  size_t runtime_objects_prealloc_size = 0;
  void *runtime_objects_prealloc_buffer = nullptr;
  if (config_.arch == Arch::cuda) {
    /*
    auto [temp_result_alloc, res] =
        llvm_device()->allocate_memory_unique({sizeof(uint64_t)});
    RW_ERROR_IF(
        res != RedwoodResult::success,
        "Failed to allocate memory for `runtime_get_memory_requirements`");
    void *temp_result_ptr = llvm_device()->get_memory_addr(*temp_result_alloc);

    redwood::runtime::runtime_get_memory_requirements((Ptr) temp_result_ptr, num_rand_states,
                                    1);
    runtime_objects_prealloc_size =
        size_t(fetch_result<uint64_t>(0, (uint64_t *)temp_result_ptr));
    temp_result_alloc.reset();
    size_t result_buffer_size = sizeof(uint64_t) * redwood_result_buffer_entries;

    RW_TRACE("Allocating device memory {:.2f} MB",
             1.0 * (runtime_objects_prealloc_size + result_buffer_size) /
                 (1UL << 20));

    runtime_objects_prealloc_buffer = preallocate_memory(
        iroundup(runtime_objects_prealloc_size + result_buffer_size,
                 redwood_page_size),
        preallocated_runtime_objects_allocs_);

    *result_buffer_ptr =
        (uint64_t *)((uint8_t *)runtime_objects_prealloc_buffer +
                     runtime_objects_prealloc_size);
    */
  } else {
    *result_buffer_ptr = (uint64_t *)HostMemoryPool::get_instance().allocate(
        sizeof(uint64_t) * redwood_result_buffer_entries, 8);
  }

  RW_TRACE("Launching runtime_initialize");

  auto *host_memory_pool = &HostMemoryPool::get_instance();
  redwood::runtime::runtime_initialize((Ptr)*result_buffer_ptr, (Ptr) host_memory_pool,
                     runtime_objects_prealloc_size,
                     (Ptr) runtime_objects_prealloc_buffer, num_rand_states,
                     (void *)&host_allocate_aligned, (void *)std::printf,
                     (void *)std::vsnprintf);

  RW_TRACE("LLVMRuntime initialized (excluding `root`)");
  llvm_runtime_ = fetch_result<void *>(redwood_result_buffer_ret_value_id,
                                       *result_buffer_ptr);
  RW_TRACE("LLVMRuntime pointer fetched");

  // Preallocate for runtime memory and update to LLVMRuntime
  if (config_.arch == Arch::cuda || config_.arch == Arch::amdgpu) {
    /*
    if (!use_device_memory_pool()) {
      preallocate_runtime_memory();
    }
    */
  }

  if (config_.arch == Arch::cuda) {
    RW_TRACE("Initializing {} random states using CUDA", num_rand_states);
    /*
    redwood::runtime::runtime_initialize_rand_states_cuda(config_.saturating_grid_dim,
        config_.max_block_dim, 0, llvm_runtime_, starting_rand_state);
      */  
  } else {
    RW_TRACE("Initializing {} random states (serially)", num_rand_states);
    auto llvm_runtime = static_cast<redwood::runtime::LLVMRuntime*>(llvm_runtime_);
    redwood::runtime::runtime_initialize_rand_states_serial(llvm_runtime, starting_rand_state);
  }

  if (arch_use_host_memory(config_.arch)) {
    auto llvm_runtime = static_cast<redwood::runtime::LLVMRuntime*>(llvm_runtime_);
    redwood::runtime::LLVMRuntime_initialize_thread_pool(llvm_runtime, thread_pool_.get(),
                                       (void *)ThreadPool::static_run);
    /*
    runtime_jit->call<void *, void *>("LLVMRuntime_set_assert_failed",
                                      llvm_runtime_,
                                      (void *)assert_failed_host);
                                      */
  }
  /*
  if (arch_is_cpu(config_.arch) && (profiler != nullptr)) {
    // Profiler functions can only be called on CPU kernels
    runtime_jit->call<void *, void *>("LLVMRuntime_set_profiler", llvm_runtime_,
                                      profiler);
    runtime_jit->call<void *, void *>(
        "LLVMRuntime_set_profiler_start", llvm_runtime_,
        (void *)&KernelProfilerBase::profiler_start);
    runtime_jit->call<void *, void *>(
        "LLVMRuntime_set_profiler_stop", llvm_runtime_,
        (void *)&KernelProfilerBase::profiler_stop);
  }
  */
}

SNodeTree* LlvmRuntimeExecutor::add_snode_tree(std::unique_ptr<SNode> root){
  const int id = allocate_snode_tree_id();
  auto tree = std::make_unique<SNodeTree>(id, std::move(root));
  tree->root()->set_snode_tree_id(id);
  //cache_data_->fields.at(id);
  /*
  if (compile_only) {
    program_impl_->compile_snode_tree_types(tree.get());
  } else {
    */
    materialize_snode_tree(tree.get(), result_buffer_);
  /*
  }
  */
  if (id < snode_trees_.size()) {
    snode_trees_[id] = std::move(tree);
  } else {
    RW_ASSERT(id == snode_trees_.size());
    snode_trees_.push_back(std::move(tree));
  }
  return snode_trees_[id].get();
}
void LlvmRuntimeExecutor::materialize_snode_tree(SNodeTree *tree,
                                      uint64_t *result_buffer_ptr){
                    
    int snode_tree_id = tree->id();
    int root_id = tree->root()->id;

    std::unique_ptr<LlvmStructCompiler> struct_compiler = std::unique_ptr<LlvmStructCompiler>();
    struct_compiler->collect_snodes(*tree->root());
    cache_field(snode_tree_id, root_id, *struct_compiler);
    initialize_llvm_runtime_snodes(cache_data_->fields.at(snode_tree_id), result_buffer_ptr);
}
void LlvmRuntimeExecutor::destroy_snode_tree(SNodeTree *snode_tree) {
  //get_llvm_context()->delete_snode_tree(snode_tree->id());
    if (cache_data_->fields.find(snode_tree->id()) != cache_data_->fields.end())
      cache_data_->fields.erase(snode_tree->id());
  snode_tree_buffer_manager_->destroy(snode_tree);
}

void LlvmRuntimeExecutor::cache_field(int snode_tree_id, int root_id, const LlvmStructCompiler& struct_compiler){
    if (cache_data_->fields.find(snode_tree_id) != cache_data_->fields.end()) {
    // [TODO] check and update the Cache, instead of simply return.
    return;
  }

  LlvmOfflineCache::FieldCacheData ret;
  ret.tree_id = snode_tree_id;
  ret.root_id = root_id;
  ret.root_size = struct_compiler.root_size;

  const auto &snodes = struct_compiler.snodes;
  for (size_t i = 0; i < snodes.size(); i++) {
    LlvmOfflineCache::FieldCacheData::SNodeCacheData snode_cache_data;
    snode_cache_data.id = snodes[i]->id;
    snode_cache_data.type = snodes[i]->type;
    snode_cache_data.cell_size_bytes = snodes[i]->cell_size_bytes;
    snode_cache_data.chunk_size = snodes[i]->chunk_size;

    ret.snode_metas.emplace_back(std::move(snode_cache_data));
  }

  cache_data_->fields[snode_tree_id] = std::move(ret);
}
int LlvmRuntimeExecutor::allocate_snode_tree_id() {
  if (free_snode_tree_ids_.empty()) {
    return snode_trees_.size();
  } else {
    int id = free_snode_tree_ids_.top();
    free_snode_tree_ids_.pop();
    return id;
  }
}

Device *LlvmRuntimeExecutor::get_compute_device() {
  return device_.get();
}

LLVMRuntime *LlvmRuntimeExecutor::get_llvm_runtime() {
  return static_cast<LLVMRuntime *>(llvm_runtime_);
}
/*
void LlvmRuntimeExecutor::init_runtime_jit_module(
    std::unique_ptr<llvm::Module> module) {
  llvm_context_->init_runtime_module(module.get());
  runtime_jit_module_ = create_jit_module(std::move(module));
}
*/

}  // namespace redwood::lang
