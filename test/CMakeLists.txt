cmake_minimum_required(VERSION 3.10)
set(This reorderingTests)

set(Sources
    branchAndBoundTests.cpp
)

add_executable(${This} ${Sources})
target_link_libraries(${This} PUBLIC 
    gtest_main
    test_lib
)
enable_testing()
add_test(
    NAME ${This}
    COMMAND ${This}
)