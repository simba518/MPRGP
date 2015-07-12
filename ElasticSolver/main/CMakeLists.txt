SET(APP_NAME elasticsolver_test)

AUX_SOURCE_DIRECTORY(. elasticsolver_test_src)

ADD_EXECUTABLE(${APP_NAME} 
  ${elasticsolver_test_src}
  )

TARGET_LINK_LIBRARIES(${APP_NAME}
  
  ${PROJECT_SOURCE_DIR}/lib/libelastic_solver.a
  /home/simba/Workspace/CommonFile/build/libDEC.so
  ${Boost_LIBRARIES}
  )