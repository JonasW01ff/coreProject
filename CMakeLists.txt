cmake_minimum_required(VERSION 3.26)
project(coreProject)

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES main.cpp)

if (APPLE)
    add_subdirectory(metalcpp)  # Library definition
endif()

include_directories(src)
add_subdirectory(src)

add_subdirectory(test)
