SET(APP_NAME collision_handle)

AUX_SOURCE_DIRECTORY(. collision_handle_src)

ADD_EXECUTABLE(${APP_NAME}
  ${collision_handle_src}
  )

TARGET_LINK_LIBRARIES(${APP_NAME}

  /home/simba/mosek/7/tools/platform/linux64x86/bin/libmosek64.so
  /usr/local/lib/libutility_common.a

  ${ELASTIC_SOLVER_LIB}
  ${COMMON_FILE_LIB}
  ${Boost_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
  )


SET(LIBRARY_NAME fem_collision)
SET(lib_src ${collision_handle_src})
LIST(REMOVE_ITEM lib_src ./main.cpp)
ADD_LIBRARY(${LIBRARY_NAME} ${lib_src})