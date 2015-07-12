project(FEM-solver)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "-fopenmp -O3 -Wall -fPIC -std=c++0x -fpermissive ${CMAKE_CXX_FLAGS}")
# set(CMAKE_CXX_FLAGS "-fopenmp -O0 -g -ggdb -fPIC -std=c++0x -fpermissive ${CMAKE_CXX_FLAGS}")

ADD_DEFINITIONS(-DDOUBLE_PRECISION)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")
SET(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib/${CMAKE_BUILD_TYPE_DIR})
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE_DIR})

# include_directories($ENV{HOME}/usr/include/CommonFile)
include_directories(/home/simba/Workspace/CommonFile)
include_directories(/home/simba/Workspace/CommonFile/solvers)
include_directories(${PROJECT_SOURCE_DIR})
link_directories($ENV{HOME}/usr/Linux/64/gcc4/lib)

# EIGEN
find_package(Eigen3 REQUIRED)
if(EIGEN3_FOUND)
    include_directories(${EIGEN3_INCLUDE_DIR})
endif(EIGEN3_FOUND)

# BOOST
SET(Boost_USE_MULTITHREADED ON)
FIND_PACKAGE(Boost REQUIRED filesystem system)
IF(Boost_FOUND)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
ELSE(Boost_FOUND)
  MESSAGE(SEND_ERROR "Cannot find boost")
ENDIF(Boost_FOUND)

# ARPACK
find_package(Arpack REQUIRED)

# GLUT
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(GLEW REQUIRED)

file(GLOB header *.h)
file(GLOB source *.cpp)

aux_source_directory(. lib_src)
list(REMOVE_ITEM lib_src ./FEMCubatureSolver.cpp)
add_library(elastic_solver ${lib_src})

add_subdirectory(main)