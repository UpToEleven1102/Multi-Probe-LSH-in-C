cmake_minimum_required(VERSION 3.12)
project(LSH_probing C)

set(CMAKE_C_STANDARD 99)

add_executable(LSH_probing main.c)

target_link_libraries(LSH_probing m)