set(SYSTEM_REDWOOD system_redwood)
add_library(${SYSTEM_REDWOOD})
target_sources(${SYSTEM_REDWOOD}
  PRIVATE
    threading.cpp
  )

  target_compile_features(system_redwood PUBLIC cxx_std_17)

target_include_directories(${SYSTEM_REDWOOD}
  PRIVATE
    ${PROJECT_SOURCE_DIR}
  )