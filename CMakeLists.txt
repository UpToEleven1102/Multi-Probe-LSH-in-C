cmake_minimum_required(VERSION 3.12)
project(LSH_probing C)

set(CMAKE_C_STANDARD 99)

add_executable(LSH_probing main.c lib/utils.c lib/utils.h lib/lsh.c lib/lsh.h lib/lsh_probing.c lib/lsh_probing.h lib/data_structure.c lib/data_structure.h)

target_link_libraries(LSH_probing m)