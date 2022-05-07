cmake_minimum_required(VERSION 3.10)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

set(RANDOM_TEST random_test)
add_executable(${RANDOM_TEST} random_test.cpp random.cpp)
add_test(NAME ${RANDOM_TEST} COMMAND ${RANDOM_TEST})
target_link_libraries(${RANDOM_TEST} ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES} pthread)

set(vector_TEST vector_test)
add_executable(${vector_TEST} vector_test.cpp random.cpp)
add_test(NAME ${vector_TEST} COMMAND ${vector_TEST})
target_link_libraries(${vector_TEST} PUBLIC ${CMAKE_PROJECT_NAME} ${GTEST_LIBRARIES} ${GTEST_MAIN_LIBRARIES} pthread)

