#pragma once

#include "vk_mem_alloc.h"

namespace redwood::lang {
    namespace vulkan{

/**
 * @brief Global VMA allocator. Only 1 instance of VmaAllocator will be created
 * for this project.
 */
extern VmaAllocator g_allocator;
    }

}  // namespace redwood::lang