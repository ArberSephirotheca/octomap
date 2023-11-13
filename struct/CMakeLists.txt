set(STRUCT_REDWOOD struct_redwood)
add_library(${STRUCT_REDWOOD})
target_sources(${STRUCT_REDWOOD}
  PRIVATE
    snode_tree.cpp
    snode_type.cpp
    snode.cpp
    struct.cpp
  )

target_include_directories(${STRUCT_REDWOOD}
  PRIVATE
    ${PROJECT_SOURCE_DIR}
  )