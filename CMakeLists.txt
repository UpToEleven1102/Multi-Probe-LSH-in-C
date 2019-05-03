cmake_minimum_required(VERSION 3.12)
project(LSH_probing C)

set(CMAKE_C_STANDARD 99)

add_executable(LSH_probing main.c lib/LSH_L2.c lib/LSH_L2.h lib/LSH_MAIN_L2.c lib/LSH_MAIN_L2.h)

target_link_libraries(LSH_probing m)