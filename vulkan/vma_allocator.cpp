// In exactly one translation unit, define the following macro before including
// the library. This will also define the VMA global allocator object.

#define VMA_IMPLEMENTATION
#define VMA_DEDICATED_ALLOCATION 0
#define VMA_DYNAMIC_VULKAN_FUNCTIONS 0

#include "vulkan/vma_allocator.h"

#include "vk_mem_alloc.h"

namespace redwood::lang {
    namespace vulkan{
VmaAllocator g_allocator;
    }
}