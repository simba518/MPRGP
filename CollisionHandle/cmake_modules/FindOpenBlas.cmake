# -Try to find OpenBlas
#
# The following are set after configuration is done: 
#  OPENBLAS_FOUND
#  OPENBLAS_LIBRARIES

SET(OPENBLAS_POSSIBLE_LIBPATHS

  /usr/lib/openblas-base
  /usr/local/lib/openblas-base
  /usr/lib64/openblas-base
  /usr/local/lib64/openblas-base

  /usr/lib
  /usr/local/lib
  /usr/lib64
  /usr/local/lib
  )

FIND_LIBRARY(OPENBLAS_LIBRARIES

  NAMES libopenblas.so
  PATHS ${OPENBLAS_POSSIBLE_LIBPATHS}
  )

IF(NOT OPENBLAS_LIBRARIES)
  MESSAGE(WARNING "OpenBlas library not found.")
  SET(OPENBLAS_FOUND FALSE)
ENDIF(NOT OPENBLAS_LIBRARIES)

IF(OPENBLAS_LIBRARIES)
  MESSAGE(STATUS "Found OpenBlas: " ${OPENBLAS_LIBRARIES})
  SET(OPENBLAS_FOUND TRUE)
ELSE(OPENBLAS_LIBRARIES)
  SET(OPENBLAS_FOUND FALSE)
ENDIF(OPENBLAS_LIBRARIES)
# MARK_AS_ADVANCED(
#   OPENBLAS_LIBRARIES
#   OPENBLAS_FOUND
#   )