set(VULKAN_REDWOOD vulkan_redwood)
add_library(${VULKAN_REDWOOD})
target_sources(${VULKAN_REDWOOD}
  PRIVATE
    buffer.cpp
    base_engine.cpp
    device_memory_pool.cpp
    vma_allocator.cpp
  )

target_include_directories(${VULKAN_REDWOOD}
  PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/extern/spdlog/include
    ${PROJECT_SOURCE_DIR}/extern/VulkanMemoryAllocator/include
    ${PROJECT_SOURCE_DIR}/extern/vk-bootstrap/src
    ${PROJECT_SOURCE_DIR}/extern/vk-bootstrap/src
  )

  target_link_libraries(${VULKAN_REDWOOD} vk-bootstrap::vk-bootstrap)