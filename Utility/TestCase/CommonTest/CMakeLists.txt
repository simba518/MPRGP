SET(APPNAME CommonTest)

AUX_SOURCE_DIRECTORY(. CommonTest_src)

ADD_EXECUTABLE(${APPNAME} ${CommonTest_src})

TARGET_LINK_LIBRARIES( ${APPNAME} 
  ${LIBDIR}/libutility_common.a
  ${Boost_LIBRARIES}
  )
