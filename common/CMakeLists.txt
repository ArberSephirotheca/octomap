set(COMMON_REDWOOD common_redwood)
set(CMAKE_CXX_STANDARD 17)
add_library(${COMMON_REDWOOD})
target_sources(${COMMON_REDWOOD}
  PRIVATE
    host_memory_pool.cpp
    unified_allocator.cpp
    logging.cpp
    core.cpp
    device.cpp
    dynamic_loader.cpp
    arch.cpp
  )
  target_compile_features(common_redwood PUBLIC cxx_std_17)

  find_package(CUDA REQUIRED)
  include_directories(${CUDA_INCLUDE_DIRS})

target_include_directories(${COMMON_REDWOOD}
  PUBLIC
    ${PROJECT_SOURCE_DIR}/extern/spdlog/include
  PRIVATE
    ${PROJECT_SOURCE_DIR}
  )

  target_compile_definitions(${COMMON_REDWOOD} PUBLIC)

  target_include_directories(${COMMON_REDWOOD}
    PRIVATE
    )

#target_link_libraries(${COMMON_REDWOOD} PRIVATE redwood_util)