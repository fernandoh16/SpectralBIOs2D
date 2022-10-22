# Project Spectral BIOs
cmake_minimum_required(VERSION 2.8)
project(UQ_BIOs)

set(SOURCE_DIR "/Users/FernandoH/Documents/ETH/Codes/UQ_BIOs/Sources")
#set(TMPDIR "/tmp")
include_directories(${SOURCE_DIR})

# Sub-directories
add_subdirectory(Tests)
add_subdirectory(Sources)

# Global options
add_definitions(-std=c++11)
add_definitions(-O3)
add_definitions(-DNDEBUG)
add_definitions(-g)
# add_definitions(-DBOOST_RESULT_OF_USE_DECLTYPE) # workaround for GCC

############
## Packages
############

# find_package(Boost)
# message(STATUS "NTL_INCLUDE_DIR = '${Boost_INCLUDE_DIRS}'" )

## config path for module finder scripts
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${SOURCE_DIR}/cmake/Modules/")
## MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
## Boost
find_package(Boost REQUIRED COMPONENTS program_options serialization filesystem system iostreams mpi)
include_directories(${Boost_INCLUDE_DIRS})
# set(CMAKE_CXX_FLAGS "-lboost_program_options")
## Eigen
find_package(Eigen3 3.1.2)
include_directories(${EIGEN3_INCLUDE_DIR})
## config path for module finder scripts
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${gMLQMC_SOURCE_DIR}/cmake/Modules/")
## NTL
find_package(NTL REQUIRED)
include_directories(${NTL_INCLUDE_DIRS})
## ZLIB
find_package(ZLIB)

# ===== Serial =====
# Add separate target for each test
foreach(Test ${UQ_BIOs_Tests})
	#include_directories(${Boost_INCLUDE_DIRS})
    add_executable(${Test} ${CMAKE_CURRENT_LIST_DIR}/Tests/${Test}.cpp ${SOURCE_FILES} ${HEADER_FILES})
    target_include_directories(${Test} PRIVATE ${BOOST_INCLUDE_DIRS})
    target_link_libraries(${Test} fftw3)
    target_link_libraries(${Test} ${Boost_LIBRARIES})
    target_link_libraries(${Test} ${NTL_LIBRARIES})
    target_link_libraries(${Test} -lz)
endforeach()
