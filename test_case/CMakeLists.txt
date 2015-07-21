SET(APPNAME collisionhandle_test)

AUX_SOURCE_DIRECTORY(. collisionhandle_test_src)

ADD_EXECUTABLE(${APPNAME} 
  # ${collisionhandle_test_src}
  main.cpp
  # ResultsToObj.cpp
  # ConvergencyCompare.cpp
  MPRGPTest.cpp
  GeneralDecoupledCompare.cpp
  # VtkToAbq.cpp
  )

TARGET_LINK_LIBRARIES( ${APPNAME}

  ${LIBRARY_OUTPUT_PATH}/libelasticsolver.a
  ${LIBRARY_OUTPUT_PATH}/libmprgp.a
  ${Boost_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
  )
