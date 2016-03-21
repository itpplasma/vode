### Basic settings
cmake_minimum_required (VERSION 2.6)

### Define the project
project (vode)
enable_language(Fortran)

### Set compiler flags
set(CMAKE_Fortran_MODULE_DIRECTORY ./OBJS)

set(CMAKE_Fortran_FLAGS "-g -Wall -Wno-unused-dummy-argument -O3")

### Define executable
add_library(vode
  dvode_f90_m.f90
)
