SET(LIBRARY_NAME utility_ui)

AUX_SOURCE_DIRECTORY(./Shader shader_src)

SET(NONE_QT_CPP

  RecordAndPlay/EventSerialization.cpp
  Render/SelfRenderRect.cpp
  Render/FileDialog.cpp
  Render/MeshRender.cpp
  Render/VolNodeGroupRender.cpp
  Render/DragGroupSelDraw.cpp
  Render/AxisTorus.cpp
  ${shader_src}
)

SET(UI_SOURCES

  RecordAndPlay/QInputEventRecorder.cpp
  RecordAndPlay/RecordAndReplayWrapper.cpp
  RecordAndPlay/QInputEventRecorderWidget.cpp
  Manipulation/SelectCtrl.cpp
  Manipulation/DragCtrl.cpp
  Animation/AniCtrl.cpp
  Manipulation/ConTrackBall.cpp
  Manipulation/Manipulatoion.cpp
  Render/QGLViewerExt.cpp
)

SET(UI_HEADERS 

  RecordAndPlay/QInputEventRecorder.h
  RecordAndPlay/RecordAndReplayWrapper.h
  RecordAndPlay/QInputEventRecorderWidget.h
  Manipulation/Manipulatoion.h
  Manipulation/SelectCtrl.h
  Manipulation/DragCtrl.h
  Manipulation/VolObjCtrl.h
  Animation/AniCtrl.h
  Animation/AniSliderBar.h
  Render/QGLViewerExt.h
  Manipulation/ConTrackBall.h
)

SET(UI_FORMS 
  RecordAndPlay/record_replay.ui
)

QT4_WRAP_CPP(UI_HEADERS_MOC ${UI_HEADERS})

SET(CMAKE_CURRENT_BINARY_DIR_RECORDED ${CMAKE_CURRENT_BINARY_DIR})
SET(CMAKE_CURRENT_BINARY_DIR  ${CMAKE_CURRENT_LIST_DIR}/RecordAndPlay)
QT4_WRAP_UI(UI_FORMS_HEADERS ${UI_FORMS})
SET(CMAKE_CURRENT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR_RECORDED})

ADD_LIBRARY(${LIBRARY_NAME}

  ${NONE_QT_CPP}
  ${UI_SOURCES}
  ${UI_HEADERS_MOC}
  ${UI_FORMS_HEADERS}
)

# install
FILE(GLOB header ${header} ./Animation/*.h ./Manipulation/*.h ./Render/*.h ./Shader/*.h ./RecordAndPlay/*.h ./Shader/*.glsl CACHE INTERNAL)
INSTALL(FILES ${header} DESTINATION ${HEADFILE_DEST})
INSTALL(TARGETS ${LIBRARY_NAME} DESTINATION ${LIB_DEST})
