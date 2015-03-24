cmake_minimum_required(VERSION 2.6)

# Default variables before reading command line arguments
# Installation default directory
# Does not work in empty 'build' directory
#set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_LIST_DIR}")



# Project name
project(ATHENA)


# Debug info
set(CMAKE_VERBOSE_MAKEFILE true)

# Our own module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")


# CMake packages
include(FindPkgConfig)


# Libraries

# Optional package cfitsio
find_package(CFITSIO)
if(CFITSIO_FOUND)
	include_directories(${CFITSIO_INCLUDE_DIR})
	set(LIBS ${LIBS} ${CFITSIO_LIBRARIES})
    	add_definitions(-D_WITH_FITS=1) 		# To be use in source files
endif()

# gsl is required only for venice
pkg_check_modules(GSL REQUIRED gsl)
link_directories(/opt/local/lib)


# Compiler flags
set(CMAKE_C_FLAGS "-Wall -Wuninitialized -pedantic -DNO_FILE_MSG -std=c99 -O3")

# Libraries
set(LIBS ${LIBS} m)


if(NOT DEFINED CMAKE_CURRENT_LIST_DIR)
    # cmake v < 2.8.3
    set(CMAKE_CURRENT_LIST_DIR ${CMAKE_CURRENT_SOURCE_DIR})
endif(NOT DEFINED CMAKE_CURRENT_LIST_DIR)
set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_LIST_DIR}")
message("Installation directory set to ${CMAKE_INSTALL_PREFIX}")


# Subdirectories
add_subdirectory(src)


# Tests
enable_testing()
message("Creating test module \"test_suite\" as " ${CMAKE_INSTALL_PREFIX} "/bin/test_suite_athena.py")
# Run tests 1 | 2 | 4 | 8 = 15
add_test(test_suite ${CMAKE_INSTALL_PREFIX}/bin/test_suite_athena.py -t 15 -P ${CMAKE_INSTALL_PREFIX}/bin -T ${CMAKE_INSTALL_PREFIX}/test)

