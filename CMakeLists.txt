cmake_minimum_required(VERSION 3.26)
project(coreProject)

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES main.cpp)

add_subdirectory(metalcpp)  # Library definition

include_directories(src)
add_subdirectory(src)


add_subdirectory(test)