SET(APPNAME tempt)

AUX_SOURCE_DIRECTORY(. tempt_src)

ADD_EXECUTABLE(${APPNAME} ${tempt_src})

TARGET_LINK_LIBRARIES( ${APPNAME} 

  ${LIBDIR}/libutility_math.a
  ${LIBDIR}/libutility_elastic.a
  ${LIBDIR}/libutility_common.a
  /usr/local/lib/libalglib.a
  ${ARPACK_LIBRARIES}
  ${SUPERLU2.0_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${CASADI_LIBRARIES}
  ${IPOPT_LIBRARIES}
  ${Boost_LIBRARIES}
  ${LAPACK_LIBRARIES}
  libpthread.so
  # /usr/lib/i386-linux-gnu/libdl.so
  )
