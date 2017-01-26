# CMake project file for fitter

cmake_minimum_required (VERSION 2.6)
project (fitter)

set (CMAKE_Fortran_COMPILER "ifort")
#set (CMAKE_C_COMPILER "icc")
#set (CMAKE_CXX_COMPILER "icpc")

enable_language (Fortran)
#enable_language (C)

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

#set (Fortran_COMPILER_NAME "ifort")

# default installation
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the NORMA directory."
      FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
message ("Fortran compiler: " ${Fortran_COMPILER_NAME})

if (Fortran_COMPILER_NAME STREQUAL "gfortran")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")

elseif (Fortran_COMPILER_NAME STREQUAL "ifort")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -fpp -fopenmp")
  #set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
elseif (Fortran_COMPILER_NAME STREQUAL "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME STREQUAL "gfortran")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME STREQUAL "gfortran")

set (SOURCES 
  "/home/mbog002/vspheroid-abm/src/vspheroid-fitter.f90" 
  "/home/mbog002/vspheroid-abm/src/precmod.f90" 
  "/home/mbog002/vspheroid-abm/src/stringmod.f90" 
)

include_directories ("../Build")

set (VSPHEROID_DIR "/projects/uoa00014/vspheroid-abm/DLL2.0/")
set (ITSOL_GLUE_DIR "/home/mbog002/itsol_glue/")
set (ITSOL_DIR "/home/mbog002/ITSOL_2/LIB/")
# build executables

add_executable (fitter ${SOURCES})
target_link_libraries(fitter ${VSPHEROID_DIR}libvspheroid-abm.so ${ITSOL_GLUE_DIR}libitsol_glue.a ${ITSOL_DIR}libitsol.a libblas.so liblapack.so)
set_target_properties(fitter PROPERTIES LINKER_LANGUAGE Fortran)

