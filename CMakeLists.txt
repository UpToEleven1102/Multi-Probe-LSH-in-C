cmake_minimum_required(VERSION 3.12)
project(LSH_probing C)

set(CMAKE_C_STANDARD 99)

add_executable(LSH_probing main.c lib/utils.c lib/utils.h lib/lsh.c lib/lsh.h)

target_link_libraries(LSH_probing m)