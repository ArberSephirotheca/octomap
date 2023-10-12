
set(CUDA_REDWOOD cuda_redwood)
add_library(${CUDA_REDWOOD})
target_sources(${CUDA_REDWOOD}
  PRIVATE
    cuda_device.cpp
    cuda_context.cpp
    cuda_driver.cpp
  )

find_package(CUDA REQUIRED)
include_directories(${CUDA_INCLUDE_DIRS})

target_include_directories(${CUDA_REDWOOD}
  PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/external/spdlog/include
  )

target_link_libraries(${CUDA_REDWOOD} PRIVATE cuda) #interop_rhi)