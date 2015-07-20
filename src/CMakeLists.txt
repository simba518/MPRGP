SET(APP_NAME collision_handle)
AUX_SOURCE_DIRECTORY(./Utility utility_src)
AUX_SOURCE_DIRECTORY(./CommonFile commonfile_src)
AUX_SOURCE_DIRECTORY(./CollisionHandle/ collisionhandle_src)
AUX_SOURCE_DIRECTORY(./MPRGPSolver/SimpleQP qp_src)

ADD_EXECUTABLE(${APP_NAME}
  ./CollisionHandle/main.cpp
  ${utility_src}
  ${commonfile_src}
  ${collisionhandle_src}
  ${qp_src}
 )

TARGET_LINK_LIBRARIES(${APP_NAME}

  ${LIBRARY_OUTPUT_PATH}/libelasticsolver.a
  ${Boost_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
  )

SET(LIBRARY_NAME mprgp)
# SET(lib_src ${utility_src} ${commonfile_src} ${collisionhandle_src} ${qp_src})
SET(lib_src ${qp_src})
LIST(REMOVE_ITEM lib_src ./CollisionHandle/main.cpp)
ADD_LIBRARY(${LIBRARY_NAME} ${lib_src})

ADD_SUBDIRECTORY(ElasticSolver)
