# ./runtime/llvm/CMakeLists.txt

add_library(llvm_runtime)
target_sources(llvm_runtime
  PRIVATE
    locked_task.cpp
    node_dynamic.cpp
    node_pointer.cpp
    runtime.cpp
    llvm_runtime_executor.cpp
    #llvm_offline_cache.cpp
    #llvm_context.cpp
    snode_tree_buffer_manager.cpp
    llvm_struct_compiler.cpp
    #kernel_launcher.cpp
  )

  target_compile_features(llvm_runtime PUBLIC cxx_std_17)
  find_package(Eigen3 REQUIRED)
target_include_directories(llvm_runtime
  PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/external/spdlog/include
    #${LLVM_INCLUDE_DIRS}
  )
  find_package(CUDA REQUIRED)
include_directories(${CUDA_INCLUDE_DIRS})

target_link_libraries(llvm_runtime PRIVATE cuda_redwood system_redwood program_redwood common_redwood)


#target_link_libraries(llvm_runtime PRIVATE ${llvm_ptx_libs})


#target_link_libraries(llvm_runtime PRIVATE ti_device_api)