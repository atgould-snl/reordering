cmake_minimum_required(VERSION 3.10)
set(This reordering)
project(${this} LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Specify the path to the Trilinos installation
set(Trilinos_DIR /fgs/atgould/teko_trilinos/trilinos-builds/installs/ifpack2/lib64/cmake/Trilinos/)

# Find Trilinos
find_package(Trilinos REQUIRED)

# Import Trilinos targets
include_directories(${Trilinos_INCLUDE_DIRS})
link_directories(${Trilinos_LIBRARY_DIRS})

set(Sources
    exhaustive.cpp
    branchAndBound.cpp
)
add_executable(${This} main.cpp ${Sources})
# Link against Trilinos and Kokkos
target_link_libraries(${This} ${Trilinos_LIBRARIES})

# Testing
include(CTest)
enable_testing()
add_subdirectory(googletest)
add_subdirectory(test)
add_library(test_lib STATIC main.cpp ${Sources})
target_link_libraries(test_lib ${Trilinos_LIBRARIES})


