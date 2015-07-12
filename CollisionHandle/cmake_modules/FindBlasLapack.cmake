# -Try to find both BLAS and LAPACK libraries.
#
# This file is similar to FindLAPACK.cmake that providging by the CMake. In
# particular, it extend the FindLAPACK.cmake library to support more optimal
# blas libraries, including OpenBlas and GotoBlas. More details could be found
# at FindLAPACK.cmake, and FindLapack.cmake. And to find more types of BLAS,
# please refer to FindBLAS.cmake.
# 
# Example:
# SET(BLA_VENDOR "OpenBlas")
# FIND_PACKAGE(BlasLapack REQUIRED)
# 
# Valid values for BLA_VENDOR including: 
# Intel(mkl), ACML,Apple, NAS, Generic, OpenBlas, Goto.
# 
# The following are set after configuration is done: 
#  LAPACK_FOUND
#  LAPACK_LIBRARIES

IF ((${BLA_VENDOR} MATCHES "OpenBlas") OR (${BLA_VENDOR} MATCHES "Goto"))

  IF(${BLA_VENDOR} MATCHES "Goto")
	FIND_PACKAGE(GotoBlas )
	IF(GOTOBLAS_FOUND)
	  SET(LAPACK_FOUND ${GOTOBLAS_FOUND})
	ENDIF(GOTOBLAS_FOUND)

  ELSE  (${BLA_VENDOR} MATCHES "Goto")
	FIND_PACKAGE(OpenBlas)
	SET( LAPACK_FOUND ${OPENBLAS_FOUND} )
  ENDIF (${BLA_VENDOR} MATCHES "Goto")

  IF ( LAPACK_FOUND )
	FIND_PACKAGE(Lapack)
	SET( LAPACK_FOUND Lapack_FOUND )
  ENDIF ( LAPACK_FOUND )

  IF ( LAPACK_FOUND )
	IF (${BLA_VENDOR} MATCHES "Goto")
	  SET ( LAPACK_LIBRARIES ${Lapack_LIBRARIES} ${GOTOBLAS_LIBRARIES} )
	ELSE (${BLA_VENDOR} MATCHES "Goto")
	  SET ( LAPACK_LIBRARIES ${Lapack_LIBRARIES} ${OPENBLAS_LIBRARIES} )
	ENDIF(${BLA_VENDOR} MATCHES "Goto")
  ENDIF ( LAPACK_FOUND )

ELSE ((${BLA_VENDOR} MATCHES "OPENBLAS") OR (${BLA_VENDOR} MATCHES "Goto"))
  FIND_PACKAGE(BLAS REQUIRED)
  FIND_PACKAGE(LAPACK REQUIRED)
  MESSAGE(STATUS "Found Generic lapack: " ${LAPACK_LIBRARIES})
  MESSAGE(STATUS "Found Generic blas: " ${BLAS_LIBRARIES})
ENDIF((${BLA_VENDOR} MATCHES "OpenBlas") OR (${BLA_VENDOR} MATCHES "Goto"))
