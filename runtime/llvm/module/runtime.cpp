#pragma once

#include <atomic>
#include <cstdint>
#include <cmath>
#include <cstdarg>
#include <cstdlib>
#include <algorithm>
#include <type_traits>
#include <cstring>

#include "inc/constants.h"
#include "common/arithmetic.h"
#include "util/bit.h"
#include "inc/cuda_kernel_utils.inc.h"
//#include "llvm-12/llvm/IR/Instructions.h"

namespace redwood{
  namespace runtime{

struct PhysicalCoordinates {
  int32_t val[redwood_max_num_indices];
};

struct RuntimeContext;


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
                

struct NodeManager;
struct LLVMRuntime;
struct ListManager;


// Common Attributes
struct StructMeta {
  int32_t snode_id;
  std::size_t element_size;
  int64_t max_num_elements;

  Ptr (*lookup_element)(Ptr, Ptr, int i);

  Ptr (*from_parent_element)(Ptr);

  char (*is_active)(Ptr, Ptr, int i);

  int32_t (*get_num_elements)(Ptr, Ptr);

  void (*refine_coordinates)(PhysicalCoordinates *inp_coord,
                             PhysicalCoordinates *refined_coord,
                             int index);

  RuntimeContext *context;
};

extern "C" {

struct Element {
  Ptr element;
  int loop_bounds[2];
  PhysicalCoordinates pcoord;
};

struct RandState {
  uint32_t x;
  uint32_t y;
  uint32_t z;
  uint32_t w;
  int32_t lock;
};

void initialize_rand_state(RandState *state, uint32_t i) {
  state->x = 123456789 * i * 1000000007;
  state->y = 362436069;
  state->z = 521288629;
  state->w = 88675123;
  state->lock = 0;
}
}


// TODO: there are many int32_t types in this class, which may be an issue if there
// are >= 2 ** 31 elements.
struct ListManager {
  static constexpr std::size_t max_num_chunks = 128 * 1024;
  Ptr chunks[max_num_chunks];
  std::size_t element_size{0};
  std::size_t max_num_elements_per_chunk;
  int32_t log2chunk_num_elements;
  int32_t lock;
  int32_t num_elements;
  LLVMRuntime *runtime;

  ListManager(LLVMRuntime *runtime,
              std::size_t element_size,
              std::size_t num_elements_per_chunk)
      : element_size(element_size),
        max_num_elements_per_chunk(num_elements_per_chunk),
        runtime(runtime) {
    /*
    redwood_assert_runtime(runtime, is_power_of_two(max_num_elements_per_chunk),
                          "max_num_elements_per_chunk must be POT.");
    */
    lock = 0;
    num_elements = 0;
    log2chunk_num_elements = redwood::log2int(num_elements_per_chunk);
  }

  void append(void *data_ptr);

  int32_t reserve_new_element() {
    // TODO: conisder using llvm::AtomicRMWInst for thread safety
    auto i = num_elements+1;//llvm::AtomicRMWInst(&num_elements, 1);
    auto chunk_id = i >> log2chunk_num_elements;
    touch_chunk(chunk_id);
    return i;
  }

  template <typename T>
  void push_back(const T &t) {
    this->append((void *)&t);
  }

  Ptr allocate();

  void touch_chunk(int chunk_id);

  int32_t get_num_active_chunks() {
    int32_t counter = 0;
    for (int i = 0; i < max_num_chunks; i++) {
      counter += (chunks[i] != nullptr);
    }
    return counter;
  }

  void clear() {
    num_elements = 0;
  }

  void resize(int32_t n) {
    num_elements = n;
  }

  Ptr get_element_ptr(int32_t i) {
    return chunks[i >> log2chunk_num_elements] +
           element_size * (i & ((1 << log2chunk_num_elements) - 1));
  }

  template <typename T>
  T &get(int32_t i) {
    return *(T *)get_element_ptr(i);
  }

  Ptr touch_and_get(int32_t i) {
    touch_chunk(i >> log2chunk_num_elements);
    return get_element_ptr(i);
  }

  int32_t size() {
    return num_elements;
  }

  int32_t ptr2index(Ptr ptr) {
    auto chunk_size = max_num_elements_per_chunk * element_size;
    for (int i = 0; i < max_num_chunks; i++) {
      //redwood_assert_runtime(runtime, chunks[i] != nullptr, "ptr not found.");
      if (chunks[i] <= ptr && ptr < chunks[i] + chunk_size) {
        return (i << log2chunk_num_elements) +
               int32_t((ptr - chunks[i]) / element_size);
      }
    }
    return -1;
  }
};

struct PreallocatedMemoryChunk {
  Ptr preallocated_head = nullptr;
  Ptr preallocated_tail = nullptr;
  std::size_t preallocated_size = 0;
};

struct LLVMRuntime {
  PreallocatedMemoryChunk runtime_objects_chunk;
  PreallocatedMemoryChunk runtime_memory_chunk;

  host_allocator_type host_allocator;
  assert_failed_type assert_failed;
  host_printf_type host_printf;
  host_vsnprintf_type host_vsnprintf;
  Ptr memory_pool;

  Ptr roots[kMaxNumSnodeTreesLlvm];
  size_t root_mem_sizes[kMaxNumSnodeTreesLlvm];

  Ptr thread_pool;
  parallel_for_type parallel_for;
  ListManager *element_lists[redwood_max_num_snodes];
  NodeManager *node_allocators[redwood_max_num_snodes];
  Ptr ambient_elements[redwood_max_num_snodes];
  Ptr temporaries;
  RandState *rand_states;

  // Cross backend (CPU, CUDA, AMDGPU) runtime memory allocation
  Ptr allocate_aligned(PreallocatedMemoryChunk &memory_chunk,
                       std::size_t size,
                       std::size_t alignment,
                       bool request = false);

  // Allocate from preallocated memory (CUDA, AMDGPU)
  Ptr allocate_from_reserved_memory(PreallocatedMemoryChunk &memory_chunk,
                                    std::size_t size,
                                    std::size_t alignment);
  Ptr profiler;
  void (*profiler_start)(Ptr, Ptr);
  void (*profiler_stop)(Ptr);

  char error_message_template[redwood_error_message_max_length];
  uint64_t error_message_arguments[redwood_error_message_max_num_arguments];
  int32_t error_message_lock = 0;
  int64_t error_code = 0;

  Ptr result_buffer;
  int32_t allocator_lock;

  int32_t num_rand_states;

  int64_t total_requested_memory;

  template <typename T>
  void set_result(std::size_t i, T t) {
    static_assert(sizeof(T) <= sizeof(uint64_t));
    ((uint64_t *)result_buffer)[i] =
        redwood_union_cast_with_different_sizes<uint64_t>(t);
  }

  template <typename T, typename... Args>
  T *create(Args &&...args) {
    auto ptr = (T *)allocate_aligned(runtime_memory_chunk, sizeof(T), 4096,
                                     true /*request*/);
    new (ptr) T(std::forward<Args>(args)...);
    return ptr;
  }
};


// NodeManager of node S (hash, pointer) managers the memory allocation of S_ch
// It makes use of three ListManagers.
struct NodeManager {
  LLVMRuntime *runtime;
  int32_t lock;

  int32_t element_size;
  int32_t chunk_num_elements;
  int32_t free_list_used;

  ListManager *free_list, *recycled_list, *data_list;
  int32_t recycle_list_size_backup;

  using list_data_type = int32_t;

  NodeManager(LLVMRuntime *runtime,
              int32_t element_size,
              int32_t chunk_num_elements = -1)
      : runtime(runtime), element_size(element_size) {
    // 128K elements per chunk, by default
    if (chunk_num_elements == -1) {
      chunk_num_elements = 128 * 1024;
    }
    // Maximum chunk size = 128 MB
    while (chunk_num_elements > 1 &&
           (uint64_t)chunk_num_elements * element_size > 128UL * 1024 * 1024) {
      chunk_num_elements /= 2;
    }
    this->chunk_num_elements = chunk_num_elements;
    free_list_used = 0;
    free_list = runtime->create<ListManager>(runtime, sizeof(list_data_type),
                                             chunk_num_elements);
    recycled_list = runtime->create<ListManager>(
        runtime, sizeof(list_data_type), chunk_num_elements);
    data_list =
        runtime->create<ListManager>(runtime, element_size, chunk_num_elements);
  }

  Ptr allocate() {
    int old_cursor = (&free_list_used, 1);
    int32_t l;
    if (old_cursor >= free_list->size()) {
      // running out of free list. allocate new.
      l = data_list->reserve_new_element();
    } else {
      // reuse
      l = free_list->get<list_data_type>(old_cursor);
    }
    return data_list->get_element_ptr(l);
  }

  int32_t locate(Ptr ptr) {
    return data_list->ptr2index(ptr);
  }

  void recycle(Ptr ptr) {
    auto index = locate(ptr);
    recycled_list->append(&index);
  }

  void gc_serial() {
    // compact free list
    for (int i = free_list_used; i < free_list->size(); i++) {
      free_list->get<list_data_type>(i - free_list_used) =
          free_list->get<list_data_type>(i);
    }
    const int32_t num_unused = std::max(free_list->size() - free_list_used, 0);
    free_list_used = 0;
    free_list->resize(num_unused);

    // zero-fill recycled and push to free list
    for (int i = 0; i < recycled_list->size(); i++) {
      auto idx = recycled_list->get<list_data_type>(i);
      auto ptr = data_list->get_element_ptr(idx);
      std::memset(ptr, 0, element_size);
      free_list->push_back(idx);
    }
    recycled_list->clear();
  }
};


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
  int i = block_dim() * block_idx() + thread_idx();
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
}
}