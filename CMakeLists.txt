### Basic settings
cmake_minimum_required (VERSION 3.13)

### Define the project
project (vode)
enable_language(Fortran)

### Find OpenMP
find_package(OpenMP REQUIRED)

### Define the library
add_library(vode STATIC
  dvode_f90_m.f90
)

# Link OpenMP to the library for threadprivate directives
target_link_libraries(vode PUBLIC OpenMP::OpenMP_Fortran)

### Testing
enable_testing()

# Thread-safety test
add_executable(test_thread_safety tests/test_thread_safety.f90)
target_link_libraries(test_thread_safety vode OpenMP::OpenMP_Fortran)
add_test(NAME thread_safety COMMAND test_thread_safety)
