cmake_minimum_required(VERSION 3.10)
project(MyApplication)

# Specify the path to the Trilinos installation
set(Trilinos_DIR /fgs/atgould/teko_trilinos/trilinos-builds/installs/ifpack2/lib64/cmake/Trilinos/)

# Find Trilinos
find_package(Trilinos REQUIRED)

# Import Trilinos targets
include_directories(${Trilinos_INCLUDE_DIRS})
link_directories(${Trilinos_LIBRARY_DIRS})

add_executable(my_application main.cpp)

# Link against Trilinos and Kokkos
target_link_libraries(my_application ${Trilinos_LIBRARIES})

# Comment change