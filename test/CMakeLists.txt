cmake_minimum_required(VERSION 3.10)
set(This reorderingTests)

# Specify the path to the Trilinos installation
set(Trilinos_DIR /fgs/atgould/teko_trilinos/trilinos-builds/installs/ifpack2/lib64/cmake/Trilinos/)

# Find Trilinos
find_package(Trilinos REQUIRED)

# Import Trilinos targets
include_directories(${Trilinos_INCLUDE_DIRS})
link_directories(${Trilinos_LIBRARY_DIRS})

# Testing
include(CTest)
enable_testing()
add_subdirectory(googletest)

set(Sources
    branchAndBoundTests.cpp
)

add_executable(${This} ${Sources})
# Link against Trilinos and Kokkos
target_link_libraries(${This} gtest_main ${Trilinos_LIBRARIES})

add_test(
    NAME ${This}
    COMMAND ${This}
)