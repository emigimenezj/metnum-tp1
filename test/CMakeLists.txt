cmake_minimum_required(VERSION 3.23)
project(test)

add_subdirectory(google-test)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

# FILES SETS
file(GLOB TESTS_TP *.cpp)

add_executable(Test_Suit ${HEADERS_TP} ${SOURCES_TP} ${TESTS_TP})
target_link_libraries(Test_Suit gtest gtest_main)