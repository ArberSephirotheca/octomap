#include "runtime/llvm/node_pointer.h"
#include "runtime/llvm/runtime.h"
#include "locked_task.h"

namespace redwood::lang
{
bool PointerNode::is_representative(uint32_t mask, uint64_t value) {
    //RW_INFO("PointerNode::is_representative");
/*
#if defined(ARCH_cuda)
  // If many threads in the mask share the same value, simply
  // elect one thread to return true and let others return false.
  if (cuda_compute_capability() < 70) {
    // <= Pascal
    bool has_following_eqiv = false;
    for (int s = 1; s < 32; s++) {
      auto cond = warp_idx() + s < 32 && ((mask >> (warp_idx() + s)) & 1);
#define TEST_PEER(x) ((x) == cuda_shfl_down_sync_i32(mask, (x), s, 31))
      auto equiv = cond && TEST_PEER(i32(i64(value))) &&
                   TEST_PEER(i32((u64)value >> 32));
#undef TEST_PEER
      has_following_eqiv = has_following_eqiv || equiv;
    }
    return !has_following_eqiv;
  } else {
    // >= Volta > Pascal
    i32 equiv_mask = cuda_match_any_sync_i64(mask, i64(value));
    auto leader = cttz_i32(equiv_mask);
    return warp_idx() == leader;
  }
#else
*/
  return true;
//#endif
}

void PointerNode::Pointer_activate(int i) {
  // TODO: something wrong with the base, check it later
  Ptr base = this->runtime->element_lists[this->id+this->get_snode_tree_id()]->get_element_ptr(0);
  RW_ASSERT(base != nullptr);
  auto num_elements = Pointer_get_num_elements();
  RW_INFO("Pointer_activate: num of max elements: {}", num_elements);
  volatile Ptr lock = (Ptr)base + 8 * i;
  *lock = std::atomic<uint64_t>(0);
  //RW_INFO("Pointer_activate: lock: {}", *lock);
  RW_INFO("Pointer_activate: lock");
  volatile Ptr *data_ptr = (Ptr *)((Ptr)base + 8 * (num_elements + i));
  RW_ASSERT(data_ptr != nullptr);
  
  if (*data_ptr == nullptr) {
    RW_INFO("Pointer_activate: data_ptr == nullptr");
    // The cuda_ calls will return 0 or do noop on CPUs
    //uint32_t mask = cuda_active_mask();
    if (is_representative(0, (uint64_t)lock)) {
        
      locked_task(
          lock,
          [&] {
            auto alloc = runtime->node_allocators[this->id];
            auto allocated = (uint64_t)alloc->allocate();
            RW_INFO("Pointer_activate: allocated");
            // TODO: Not sure if we really need atomic_exchange here,
            // just to be safe.
            std::memcpy((void *)data_ptr, &allocated, sizeof(uint64_t));
           // atomic_exchange_u64((uint64_t *)data_ptr, allocated);
          },
          [&]() { return *data_ptr == nullptr; });
    }
    //warp_barrier(mask);
  }
}

void PointerNode::Pointer_deactivate(int i) {
  auto num_elements = Pointer_get_num_elements();
  Ptr lock = (Ptr) this + 8 * i;
  Ptr &data_ptr = *(Ptr *)(this + 8 * (num_elements + i));
  if (data_ptr != nullptr) {
    locked_task(lock, [&] {
      if (data_ptr != nullptr) {
        auto alloc = runtime->node_allocators[this->id];
        alloc->recycle(data_ptr);
        data_ptr = nullptr;
      }
    });
  }
}

bool PointerNode::Pointer_is_active(int i) {
  auto num_elements = Pointer_get_num_elements();
  auto data_ptr = *(Ptr *)(this + 8 * (num_elements + i));
  return data_ptr != nullptr;
}

Ptr PointerNode::Pointer_lookup_element(int i) {
  auto num_elements = Pointer_get_num_elements();
  auto data_ptr = *(Ptr *)(this + 8 * (num_elements + i));
  if (data_ptr == nullptr) {
    data_ptr = (runtime)->ambient_elements[this->id];
  }
  return data_ptr;
}
}