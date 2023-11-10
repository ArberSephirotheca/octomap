set(PROGRAM_REDWOOD program_redwood)
add_library(${PROGRAM_REDWOOD})
target_sources(${PROGRAM_REDWOOD}
  PRIVATE
    compile_config.cpp
  )

target_include_directories(${PROGRAM_REDWOOD}
  PRIVATE
    ${PROJECT_SOURCE_DIR}
  )

  