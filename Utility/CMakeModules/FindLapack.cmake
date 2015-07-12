# -Try to find LAPACK (find LAPACK only)
# 
# This file is different to the FindLAPACK.cmake that provided by the CMake, as
# the latter will return both BLAS and LAPACK library. See also
# FindBlasLapack.cmake.
# 
# The following are set after configuration is done: 
#  Lapack_FOUND
#  Lapack_LIBRARIES

SET(Lapack_POSSIBLE_LIBPATHS

  /usr/lib
  /usr/local/lib
  /usr/lib64
  /usr/local/lib
  )

FIND_LIBRARY(Lapack_LIBRARIES

  NAMES lapack
  PATHS ${Lapack_POSSIBLE_LIBPATHS}
  )

IF(NOT Lapack_LIBRARIES)
  MESSAGE(SEND_ERROR "Lapack library not found.")
ENDIF(NOT Lapack_LIBRARIES)

IF(Lapack_LIBRARIES)
  MESSAGE(STATUS "Found Lapack")
  SET(Lapack_FOUND TRUE)
ELSE(Lapack_LIBRARIES)
  SET(Lapack_FOUND FALSE)
ENDIF(Lapack_LIBRARIES)

MARK_AS_ADVANCED(

  Lapack_LIBRARIES
  Lapack_FOUND
  )