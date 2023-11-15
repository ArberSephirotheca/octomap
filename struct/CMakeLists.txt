set(STRUCT_REDWOOD struct_redwood)
add_library(${STRUCT_REDWOOD})
target_sources(${STRUCT_REDWOOD}
  PRIVATE
    snode_tree.cpp
    snode_type.cpp
    snode.cpp
    struct.cpp
  )

  target_compile_features(struct_redwood PUBLIC cxx_std_17)

target_include_directories(${STRUCT_REDWOOD}
  PRIVATE
    ${PROJECT_SOURCE_DIR}
  )
  target_link_libraries(struct_redwood PRIVATE llvm_runtime)
