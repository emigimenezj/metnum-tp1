cmake_minimum_required(VERSION 3.23)
project(c)

set(CMAKE_CXX_STANDARD 14)

# FILES SETS
file(GLOB HEADERS_TP model/*.h)
file(GLOB SOURCES_TP model/*.cpp)

add_executable(c main.cpp ${HEADERS_TP} ${SOURCES_TP})
add_subdirectory(test)
