SET(APPNAME collisionhandle_test)

AUX_SOURCE_DIRECTORY(. collisionhandle_test_src)

ADD_EXECUTABLE(${APPNAME} 
  # ${collisionhandle_test_src}
  main.cpp
  # ResultsToObj.cpp
  ConvergencyCompare.cpp
  # VtkToAbq.cpp
  )

TARGET_LINK_LIBRARIES( ${APPNAME}

  ${LIBDIR}/libfem_collision.a
  /home/simba/mosek/7/tools/platform/linux64x86/bin/libmosek64.so
  /usr/local/lib/libutility_common.a
  ${ELASTIC_SOLVER_LIB}
  ${COMMON_FILE_LIB}
  ${Boost_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
  )