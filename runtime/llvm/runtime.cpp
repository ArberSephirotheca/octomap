#include "runtime/llvm/runtime.h"
#include "struct/snode_types.h"
//#include "runtime/llvm/locked_task.h"
//#include "runtime/llvm/node_dynamic.h"
namespace redwood{
  namespace runtime{

using Ptr = uint8_t *;
using assert_failed_type = void (*)(const char *);
using host_printf_type = void (*)(const char *, ...);
using host_vsnprintf_type = int (*)(char *,
                                    std::size_t,
                                    const char *,
                                    std::va_list);
using host_allocator_type = void *(*)(void *, std::size_t, std::size_t);
using RangeForTaskFunc = void(RuntimeContext *, const char *tls, int i);
using MeshForTaskFunc = void(RuntimeContext *, const char *tls, int32_t i);
using parallel_for_type = void (*)(void *thread_pool,
                                   int splits,
                                   int num_desired_threads,
                                   void *context,
                                   void (*func)(void *, int thread_id, int i));




void initialize_struct_meta(redwood::lang::SNode* snode, StructMeta* meta){
  // TODO： unimplemented
  if(snode->type == redwood::lang::SNodeType::dynamic){
    
  }
  meta->snode_id  snd
  meta->element_size = get_type_size(snode);
}

std::size_t get_type_size(redwood::lang::SNode* snode){
  return sizeof(snode->data_type);
}
void ListManager::touch_chunk(int chunk_id) {
  if (!chunks[chunk_id]){
    auto chunk_ptr = runtime->allocate_aligned(
            runtime->runtime_memory_chunk,
            max_num_elements_per_chunk * element_size, 4096, true);
    // TODO: make it atomic
    chunks[chunk_id] = chunk_ptr;
    //std::atomic_exchange(&chunks[chunk_id], chunk_ptr);
  }
  /*
  if (!chunks[chunk_id]) {
    locked_task(&lock, [&] {
      // may have been allocated during lock contention
      if (!chunks[chunk_id]) {
        grid_memfence();
        auto chunk_ptr = runtime->allocate_aligned(
            runtime->runtime_memory_chunk,
            max_num_elements_per_chunk * element_size, 4096, true);
        atomic_exchange_u64((u64 *)&chunks[chunk_id], (u64)chunk_ptr);
      }
    });
  }
  */
}


void ListManager::append(void *data_ptr) {
  auto ptr = allocate();
  std::memcpy(ptr, data_ptr, element_size);
}

Ptr ListManager::allocate() {
  auto i = reserve_new_element();
  return get_element_ptr(i);
}


void initialize_rand_state(RandState *state, uint32_t i) {
  state->x = 123456789 * i * 1000000007;
  state->y = 362436069;
  state->z = 521288629;
  state->w = 88675123;
  state->lock = 0;
}


// [ON HOST] CPU backend
// [ON DEVICE] CUDA/AMDGPU backend
Ptr LLVMRuntime::allocate_aligned(PreallocatedMemoryChunk &memory_chunk,
                                  std::size_t size,
                                  std::size_t alignment,
                                  bool request) {
  // TODO: make it atomic
  if (request)
    total_requested_memory += size;
  //std::atomic_fetch_add(&total_requested_memory, size);
    //atomic_add_i64(&total_requested_memory, size);

  // TODO: do it later for CUDA
  /*
  if (memory_chunk.preallocated_size > 0) {
    return allocate_from_reserved_memory(memory_chunk, size, alignment);
  }
  */

  return (Ptr)host_allocator(memory_pool, size, alignment);
}
/*
#include "runtime/llvm/module/locked_task.h"
// [ONLY ON DEVICE] CUDA/AMDGPU backend
Ptr LLVMRuntime::allocate_from_reserved_memory(
    PreallocatedMemoryChunk &memory_chunk,
    std::size_t size,
    std::size_t alignment) {
  Ptr ret = nullptr;
  bool success = false;
  locked_task(&allocator_lock, [&] {
    std::size_t preallocated_head = (std::size_t)memory_chunk.preallocated_head;
    std::size_t preallocated_tail = (std::size_t)memory_chunk.preallocated_tail;

    auto alignment_bytes =
        alignment - 1 - (preallocated_head + alignment - 1) % alignment;
    size += alignment_bytes;
    if (preallocated_head + size <= preallocated_tail) {
      ret = (Ptr)(preallocated_head + alignment_bytes);
      memory_chunk.preallocated_head += size;
      success = true;
    } else {
      success = false;
    }
  });
  if (!success) {
#if ARCH_cuda
    // Here unfortunately we have to rely on a native CUDA assert failure to
    // halt the whole grid. Using a taichi_assert_runtime will not finish the
    // whole kernel execution immediately.
    __assertfail(
        "Out of CUDA pre-allocated memory.\n"
        "Consider using ti.init(device_memory_fraction=0.9) or "
        "ti.init(device_memory_GB=4) to allocate more"
        " GPU memory",
        "Redwood JIT", 0, "allocate_from_reserved_memory", 1);
#endif
  }
  //redwood_assert_runtime(this, success, "Out of pre-allocated memory");
  return ret;
}
*/
// External API
// [ON HOST] CPU backend
// [ON DEVICE] CUDA/AMDGPU backend
void runtime_memory_allocate_aligned(LLVMRuntime *runtime,
                                     std::size_t size,
                                     std::size_t alignment,
                                     uint64_t *result) {
  *result =
      redwood_union_cast_with_different_sizes<uint64_t>(runtime->allocate_aligned(
          runtime->runtime_memory_chunk, size, alignment));
}

// External API
// [ON HOST] CPU backend
// [ON DEVICE] CUDA/AMDGPU backend
void runtime_initialize(
    Ptr result_buffer,
    Ptr memory_pool,
    std::size_t
        preallocated_size,  // Non-zero means use the preallocated buffer
    Ptr preallocated_buffer,
    int32_t num_rand_states,
    void *_host_allocator,
    void *_host_printf,
    void *_host_vsnprintf) {
  // bootstrap
  auto host_allocator = (host_allocator_type)_host_allocator;
  auto host_printf = (host_printf_type)_host_printf;
  auto host_vsnprintf = (host_vsnprintf_type)_host_vsnprintf;
  LLVMRuntime *runtime = nullptr;
  Ptr preallocated_tail = preallocated_buffer + preallocated_size;
  if (preallocated_size) {
    runtime = (LLVMRuntime *)preallocated_buffer;
    preallocated_buffer +=
        redwood::iroundup(sizeof(LLVMRuntime), redwood_page_size);
  } else {
    runtime =
        (LLVMRuntime *)host_allocator(memory_pool, sizeof(LLVMRuntime), 128);
  }

  PreallocatedMemoryChunk runtime_objects_chunk;
  runtime_objects_chunk.preallocated_size = preallocated_size;
  runtime_objects_chunk.preallocated_head = preallocated_buffer;
  runtime_objects_chunk.preallocated_tail = preallocated_tail;

  runtime->runtime_objects_chunk = std::move(runtime_objects_chunk);

  runtime->result_buffer = result_buffer;
  runtime->set_result(redwood_result_buffer_ret_value_id, runtime);
  runtime->host_allocator = host_allocator;
  runtime->host_printf = host_printf;
  runtime->host_vsnprintf = host_vsnprintf;
  runtime->memory_pool = memory_pool;

  runtime->total_requested_memory = 0;

  runtime->temporaries = (Ptr)runtime->allocate_aligned(
      runtime->runtime_objects_chunk, redwood_global_tmp_buffer_size,
      redwood_page_size);

  runtime->num_rand_states = num_rand_states;
  runtime->rand_states = (RandState *)runtime->allocate_aligned(
      runtime->runtime_objects_chunk,
      sizeof(RandState) * runtime->num_rand_states, redwood_page_size);
}

void runtime_initialize_snodes(LLVMRuntime *runtime,
                               std::size_t root_size,
                               const int root_id,
                               const int num_snodes,
                               const int snode_tree_id,
                               std::size_t rounded_size,
                               Ptr ptr,
                               bool all_dense) {
  // For Metal runtime, we have to make sure that both the beginning address
  // and the size of the root buffer memory are aligned to page size.
  runtime->root_mem_sizes[snode_tree_id] = rounded_size;
  runtime->roots[snode_tree_id] = ptr;
  // runtime->request_allocate_aligned ready to use
  // initialize the root node element list
  if (all_dense) {
    return;
  }
  for (int i = root_id; i < root_id + num_snodes; i++) {
    // TODO: some SNodes do not actually need an element list.
    runtime->element_lists[i] =
        runtime->create<ListManager>(runtime, sizeof(Element), 1024 * 64);
  }
  
  Element elem;
  elem.loop_bounds[0] = 0;
  elem.loop_bounds[1] = 1;
  elem.element = runtime->roots[snode_tree_id];
  for (int i = 0; i < redwood_max_num_indices; i++) {
    elem.pcoord.val[i] = 0;
  }

  runtime->element_lists[root_id]->append(&elem);
}

// External API
// [ON HOST] CPU backend
// [ON DEVICE] CUDA/AMDGPU backend
void runtime_get_memory_requirements(Ptr result_buffer,
                                     int32_t num_rand_states,
                                     int32_t use_preallocated_buffer) {
  int64_t size = 0;

  if (use_preallocated_buffer) {
    size += redwood::iroundup(int64_t(sizeof(LLVMRuntime)), redwood_page_size);
  }

  size +=
      redwood::iroundup(int64_t(redwood_global_tmp_buffer_size), redwood_page_size);
  size += redwood::iroundup(int64_t(sizeof(RandState)) * num_rand_states,
                           redwood_page_size);

  reinterpret_cast<int64_t *>(result_buffer)[0] = size;
}

void runtime_NodeAllocator_initialize(LLVMRuntime *runtime,
                                      int snode_id,
                                      std::size_t node_size) {
  runtime->node_allocators[snode_id] =
      runtime->create<NodeManager>(runtime, node_size, 1024 * 16);
}

void runtime_allocate_ambient(LLVMRuntime *runtime,
                              int snode_id,
                              std::size_t size) {
  // Do not use NodeManager for the ambient node since it will never be garbage
  // collected.
  runtime->ambient_elements[snode_id] = runtime->allocate_aligned(
      runtime->runtime_memory_chunk, size, 128, true /*request*/);
}

void runtime_initialize_memory(LLVMRuntime *runtime,
                               std::size_t preallocated_size,
                               Ptr preallocated_buffer) {
  if (preallocated_size) {
    runtime->runtime_memory_chunk.preallocated_size = preallocated_size;
    runtime->runtime_memory_chunk.preallocated_head = preallocated_buffer;
    runtime->runtime_memory_chunk.preallocated_tail =
        preallocated_buffer + preallocated_size;
  }
}
void runtime_initialize_rand_states_cuda(LLVMRuntime *runtime,
                                         int starting_rand_state) {
  // TODO： fix the linking problem
  //int i = block_dim() * block_idx() + thread_idx();
  int i = 5;
  initialize_rand_state(&runtime->rand_states[i], starting_rand_state + i);
}
void runtime_initialize_rand_states_serial(LLVMRuntime *runtime,
                                           int starting_rand_state) {
  for (int i = 0; i < runtime->num_rand_states; i++) {
    initialize_rand_state(&runtime->rand_states[i], starting_rand_state + i);
  }
}

void LLVMRuntime_initialize_thread_pool(LLVMRuntime *runtime,
                                        void *thread_pool,
                                        void *parallel_for) {
  runtime->thread_pool = (Ptr)thread_pool;
  runtime->parallel_for = (parallel_for_type)parallel_for;
}

}
}