#include "runtime/llvm/node_dynamic.h"
#include "runtime/llvm/runtime.h"
#include "locked_task.h"

namespace redwood::lang
{

    // We need to not only update node->n, but also make sure the chunk containing
    // element i is allocated.
    void DynamicNode::Dynamic_activate(int i)
    {
        atomic_max_i32(n, i + 1);
        int chunk_start = 0;
        auto p_chunk_ptr = &ptr;
        while (true)
        {
            if (*p_chunk_ptr == nullptr)
            {
              RW_INFO("p_chunk_ptr == nullptr for id: {}", this->id);
              auto alloc = runtime->node_allocators[this->id];
              if(alloc == nullptr){
                
              }
          *p_chunk_ptr = alloc->allocate();
              /*
                locked_task(Ptr(&lock), [&]
                            {
        if (*p_chunk_ptr == nullptr) {
          auto alloc = runtime->node_allocators[this->id];
          *p_chunk_ptr = alloc->allocate();
        } });*/
            }
            if (i < chunk_start + chunk_size)
            {
                return;
            }
            p_chunk_ptr = (Ptr *)*p_chunk_ptr;
            chunk_start += chunk_size;
        }
    }
    void DynamicNode::Dynamic_deactivate()
    {
        if (n > 0)
        {
            locked_task(Ptr(&lock), [&]
                        {
      n = 0;
      auto p_chunk_ptr = &ptr;
      auto alloc = runtime->node_allocators[this->id];
      while (*p_chunk_ptr) {
        alloc->recycle(*p_chunk_ptr);
        p_chunk_ptr = (Ptr *)*p_chunk_ptr;
      }
      ptr = nullptr; });
        }
    }

    Ptr DynamicNode::Dynamic_allocate(int32_t *len)
    {
        auto i = atomic_max_i32(n, 1);
        *len = i;
        int chunk_start = 0;
        auto p_chunk_ptr = &ptr;
        while (true)
        {
            if (*p_chunk_ptr == nullptr)
            {
                locked_task(Ptr(&lock), [&]
                            {
        if (*p_chunk_ptr == nullptr) {
          auto alloc = runtime->node_allocators[this->id];
          *p_chunk_ptr = alloc->allocate();
        } });
            }
            if (i < chunk_start + chunk_size)
            {
                return *p_chunk_ptr + sizeof(Ptr) +
                       (i - chunk_start) * cell_size_bytes;
            }
            p_chunk_ptr = (Ptr *)(*p_chunk_ptr);
            chunk_start += chunk_size;
        }
        // Unreachable
        return nullptr;
    }
Ptr DynamicNode::Dynamic_lookup_element( int i) {
  if (Dynamic_is_active(i)) {
    int chunk_start = 0;
    auto chunk_ptr = ptr;
    while (true) {
      if (i < chunk_start + chunk_size) {
        auto addr =
            chunk_ptr + sizeof(Ptr) + (i - chunk_start) * cell_size_bytes;
        return addr;
      }
      chunk_ptr = *(Ptr *)chunk_ptr;
      chunk_start += chunk_size;
    }
  } else {
    return (runtime)->ambient_elements[this->id];
  }
}
  /*
  void DynamicNode::set_cell_size_bytes(std::size_t size){
    RW_INFO("DynamicNode::set_cell_size_bytes");
    cell_size_bytes = size;
    this->runtime->node_allocators[this->id]->set_element_size_bytes(size);
  }
  */
}