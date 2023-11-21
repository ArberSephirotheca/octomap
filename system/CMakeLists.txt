set(SYSTEM_REDWOOD system_redwood)
add_library(${SYSTEM_REDWOOD})
target_sources(${SYSTEM_REDWOOD}
  PRIVATE
    threading.cpp
  )

  target_compile_features(system_redwood PUBLIC cxx_std_17)

target_include_directories(${SYSTEM_REDWOOD}
  PUBLIC
  ${PROJECT_SOURCE_DIR}/extern/spdlog/include
  PRIVATE
    ${PROJECT_SOURCE_DIR}
  )