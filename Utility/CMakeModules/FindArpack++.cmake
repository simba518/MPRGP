# -Try to find Arpack++
#
#
# The following are set after configuration is done: 
# ARPACK++_FOUND
# ARPACK++_INCLUDE_DIR


FIND_PATH(ARPACK++_INCLUDE_DIR armat.h
  /usr/include/arpack++
  /usr/local/include/arpack++
)

#MESSAGE("DBG ARPACK++_INCLUDE_DIR=${ARPACK++_INCLUDE_DIR}")

IF(NOT ARPACK++_INCLUDE_DIR)
       MESSAGE(SEND_ERROR "Arpack++ include dir not found.")
ENDIF(NOT ARPACK++_INCLUDE_DIR)

IF(ARPACK++_INCLUDE_DIR)
	MESSAGE(STATUS "Using Arpack++ from ${ARPACK++_INCLUDE_DIR}")
	SET(ARPACK++_FOUND TRUE)
ELSE(ARPACK++_INCLUDE_DIR)
	SET(ARPACK++_FOUND FALSE)
ENDIF(ARPACK++_INCLUDE_DIR)

MARK_AS_ADVANCED(
  ARPACK++_INCLUDE_DIR
  ARPACK++_FOUND
)