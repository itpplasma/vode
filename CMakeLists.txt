### Basic settings
cmake_minimum_required (VERSION 3.13)

### Define the project
project (vode)
enable_language(Fortran)

### Define executable
add_library(vode STATIC
  dvode_f90_m.f90
)
