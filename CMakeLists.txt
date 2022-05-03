cmake_minimum_required(VERSION 3.12)

project(position-based-fluids CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED)

file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp)
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

add_executable(main ${headers} ${sources})
target_link_libraries(main PUBLIC Eigen3::Eigen)
