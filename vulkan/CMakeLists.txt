set(VULKAN_REDWOOD vulkan_redwood)
add_library(${VULKAN_REDWOOD})
target_sources(${VULKAN_REDWOOD}
  PRIVATE
    cpu_device.cpp
  )

target_include_directories(${VULKAN_REDWOOD}
  PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/external/spdlog/include
    ${PROJECT_SOURCE_DIR}/external/VulkanMemoryAllocator/include
  )