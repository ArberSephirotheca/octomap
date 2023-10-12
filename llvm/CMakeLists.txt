set(LLVM_REDWOOD llvm_redwood)
add_library(${LLVM_REDWOOD})
target_sources(${LLVM_REDWOOD}
  PRIVATE
    allocator.cpp
    device_memory_pool.cpp
  )

target_include_directories(${LLVM_REDWOOD}
  PRIVATE
    ${PROJECT_SOURCE_DIR}
  )

target_link_libraries(${LLVM_REDWOOD} PRIVATE common_redwood)