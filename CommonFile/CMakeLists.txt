PROJECT(CommonFile)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

# debug or release
SET(CMAKE_CXX_FLAGS "-fopenmp -O3 -Wall -fPIC -fpermissive -std=c++0x ${CMAKE_CXX_FLAGS}")
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/;${CMAKE_MODULE_PATH}")

ADD_DEFINITIONS(-DDOUBLE_PRECISION)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})
# EIGEN3
FIND_PACKAGE(Eigen3 REQUIRED)
IF(EIGEN3_FOUND)
  INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
ENDIF(EIGEN3_FOUND)

# BOOST
SET(Boost_USE_MULTITHREADED ON)
FIND_PACKAGE(Boost REQUIRED filesystem system)
IF(Boost_FOUND)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
ELSE(Boost_FOUND)
  MESSAGE(SEND_ERROR "Cannot find boost")
ENDIF(Boost_FOUND)

# add library
# INCLUDE_DIRECTORIES(${OPENCL_DIRS})

FILE(GLOB header ${PROJECT_SOURCE_DIR}/*.h  )
FILE(GLOB source ${PROJECT_SOURCE_DIR}/*.cpp)
FILE(GLOB solver_header ${PROJECT_SOURCE_DIR}/solvers/*.h  )
FILE(GLOB solver_source ${PROJECT_SOURCE_DIR}/solvers/*.cpp)

# for opencl
#SET(opencl_dirs "unknown opencl" CACHE PATH "The OpenCL include Path")
#IF(NOT DEFINED opencl_dirs)
#  MESSAGE(STATUS "WARNNING: Please define opencl_dirs!, Otherwise, cl won't be compiled!")
#ELSE(NOT DEFINED opencl_dirs)
#IF(NOT EXISTS "${opencl_dirs}")
#  MESSAGE(STATUS "WARNNING: Cannot find ${opencl_dirs}!, cl won't be compiled!")
#ENDIF(NOT EXISTS "${opencl_dirs}")
#ELSE(NOT EXISTS "${opencl_dirs}")
#  FILE(GLOB source ${source} ${PROJECT_SOURCE_DIR}/cl/*.h)
#  FILE(GLOB header ${header} ${PROJECT_SOURCE_DIR}/cl/*.cpp)
#ELSE(NOT EXISTS "${opencl_dirs}")
#ENDIF(NOT DEFINED opencl_dirs)

# build library
ADD_LIBRARY(DEC SHARED ${header} ${source} ${solver_header} ${solver_source})
TARGET_LINK_LIBRARIES(DEC ${LIBRARY_NAME} ${Boost_LIBRARIES} gomp)

INSTALL(FILES ${header}
	DESTINATION ${CMAKE_INSTALL_PREFIX}/include/CommonFile)
INSTALL(FILES ${solver_header}
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include/CommonFile/solvers)
#INSTALL(FILES ${CMAKE_SOURCE_DIR}/cl/*.h
#        DESTINATION ${CMAKE_INSTALL_PREFIX}/include/CommonFile/cl)
INSTALL(TARGETS DEC
	DESTINATION ${CMAKE_INSTALL_PREFIX}/Linux/64/gcc4/lib/)

