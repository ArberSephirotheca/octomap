cmake_minimum_required(VERSION 3.16)
project(my_app LANGUAGES CXX)

add_executable(my_app)

file(GLOB HEADER_LIST CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/include/*.hpp")

target_sources(
  my_app
  PRIVATE src/main.cpp src/binary_radix_tree.cpp
  PUBLIC ${HEADER_LIST})

target_include_directories(
  my_app
  PRIVATE ${CMAKE_SOURCE_DIR}/extern/libmorton/include
          ${CMAKE_SOURCE_DIR}/extern/cxxopts/include
  PUBLIC ${CMAKE_SOURCE_DIR}/include)

target_compile_features(my_app PUBLIC cxx_std_17)
target_compile_options(my_app PRIVATE -Wall -Wextra -Wpedantic)

find_package(realsense2 REQUIRED)

# Add Eigen
find_package(Eigen3 REQUIRED)
target_include_directories(my_app PUBLIC ${EIGEN3_INCLUDE_DIR})
message(STATUS "Eigen3 version: ${Eigen3_VERSION}")


find_package(PCL 1.10 REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

find_package(OpenMP REQUIRED)
# Add OpenMP flags to the compiler
if(OpenMP_CXX_FOUND)
    target_compile_options(my_app PRIVATE ${OpenMP_CXX_FLAGS})
endif()

# Link OpenMP library
if(OpenMP_CXX_FOUND)
    target_link_libraries(my_app PRIVATE OpenMP::OpenMP_CXX ${realsense2_LIBRARY} ${PCL_LIBRARIES})
endif()