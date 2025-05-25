cmake_minimum_required(VERSION 4.0)
project(PhysicsStuff)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# futil.h
add_library(futil STATIC lib/futil/futil.cpp)
target_include_directories(futil PUBLIC ${CMAKE_SOURCE_DIR}/lib)

# mutil
add_library(mutil STATIC
    lib/mutil/interp.cpp
    # lib/mutil/quad.cpp
)
target_include_directories(mutil PUBLIC ${CMAKE_SOURCE_DIR}/lib)

# ch03: oscillation
add_executable(oscillation ${CMAKE_SOURCE_DIR}/src/ch03/oscillation.cpp)
target_link_libraries(oscillation PRIVATE futil)
set_target_properties(oscillation PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/ch03)

# ch03: scattering
add_executable(scattering ${CMAKE_SOURCE_DIR}/src/ch03/scattering.cpp)
target_link_libraries(scattering PRIVATE futil)
set_target_properties(scattering PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/ch03)
