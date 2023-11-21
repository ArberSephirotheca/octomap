set(CPU_REDWOOD cpu_redwood)
add_library(${CPU_REDWOOD})
target_sources(${CPU_REDWOOD}
  PRIVATE
    cpu_device.cpp
  )

target_include_directories(${CPU_REDWOOD}
  PRIVATE
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/extern/eigen
    ${PROJECT_SOURCE_DIR}/extern/spdlog/include
  )