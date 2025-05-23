cmake_minimum_required(VERSION 4.0)
project(PhysicsStuff)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# util.h
add_library(util STATIC lib/util/util.cpp)
target_include_directories(util PUBLIC ${CMAKE_SOURCE_DIR}/lib)

# ch03: oscillation
add_executable(oscillation ${CMAKE_SOURCE_DIR}/src/ch03/oscillation.cpp)
target_link_libraries(oscillation PRIVATE util)
set_target_properties(oscillation PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/ch03)
