### Basic settings
cmake_minimum_required (VERSION 2.6)

### Define the project
project (vode)
enable_language(Fortran)

### Set compiler flags
#set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/OBJS)

set(CMAKE_Fortran_FLAGS "-O3")


### Define executable
add_library(vode SHARED
  dvode_f90_m.f90
)
