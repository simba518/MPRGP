SET(LIBRARY_NAME utility_cuda)
AUX_SOURCE_DIRECTORY(. CUDA_MATH_src)

###################CUDA###################
SET(BUILD_SHARED_LIBS OFF)
SET(CUDA_NVCC_FLAGS --compiler-options '-fPIC' --use_fast_math -arch=compute_13 -code=sm_13)
SET(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE ON)
CUDA_ADD_LIBRARY(${LIBRARY_NAME}
  
  ${CUDA_MATH_src}
  EleProduct.cu
  )
###################CUDA###################

# install
FILE(GLOB header ${header} *.h)
INSTALL(FILES ${header} DESTINATION ${HEADFILE_DEST})
INSTALL(TARGETS ${LIBRARY_NAME} DESTINATION ${LIB_DEST})