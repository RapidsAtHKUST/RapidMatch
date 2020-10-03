set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -g -Wall")
set(SRC_FILES nop_join.cpp hashjoins_test.cpp)
add_executable(hashjoins_test ${SRC_FILES})