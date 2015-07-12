SET(APPNAME RenderTest)

SET(NONE_QT_CPP
  main.cpp
)

SET(UI_HEADERS 
  QtEventRecorder.h
)

QT4_WRAP_CPP(UI_HEADERS_MOC ${UI_HEADERS})

ADD_EXECUTABLE(${APPNAME} 
  ${NONE_QT_CPP}
  ${UI_HEADERS_MOC}
)

TARGET_LINK_LIBRARIES( ${APPNAME} 

  ${LIBDIR}/libutility_ui.a
  ${LIBDIR}/libutility_common.a
  ${LIBDIR}/libutility_math.a
  ${QGLViewer_LIBRARIES}
  ${QT_LIBRARIES}
  ${OPENGL_LIBRARIES}
  ${Boost_LIBRARIES}
  )
