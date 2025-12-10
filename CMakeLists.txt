cmake_minimum_required(VERSION 3.13)

project(vode LANGUAGES Fortran)

find_package(OpenMP REQUIRED)

add_library(vode STATIC dvode_f90_m.f90)
target_link_libraries(vode PUBLIC OpenMP::OpenMP_Fortran)

# Testing only if not used as a subproject
# Disable tests by setting either VODE_BUILD_TESTING or BUILD_TESTING to OFF earlier
set(VODE_IS_STANDALONE OFF)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    set(VODE_IS_STANDALONE ON)
endif()

option(VODE_BUILD_TESTING "Build tests for the vode library" ${VODE_IS_STANDALONE})

if(DEFINED BUILD_TESTING)
    set(VODE_BUILD_TESTING ${BUILD_TESTING})
endif()

if(VODE_BUILD_TESTING)
    include(CTest)

    add_executable(test_thread_safety tests/test_thread_safety.f90)
    target_link_libraries(test_thread_safety vode OpenMP::OpenMP_Fortran)
    add_test(NAME thread_safety COMMAND test_thread_safety)
endif()
