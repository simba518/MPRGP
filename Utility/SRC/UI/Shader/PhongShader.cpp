#include "PhongShader.h"
#include <Log.h>
using namespace UTILITY;

pPhongShader PhongShader::p_instance;

void PhongShader::initialize(){

  if (pSubFaceGLSL){
	// pSubFaceGLSL should be initialized only once.
	return ;
  }

  const GLenum err = glewInit();
  if (err != GLEW_OK){

	ERROR_LOG("failed to initialize the glew." << glewGetErrorString(err));
  }else{

	if (!glewIsSupported("GL_VERSION_2_0 GL_VERSION_1_5 GL_ARB_multitexture GL_ARB_vertex_buffer_object")) {

	  ERROR_LOG("failed to initialize the glew, as certern API of glew is not supported");
	}else{
	  
	  const std::string SHADER_PATH = "/usr/local/include/utility/";
	  const std::string vertexShader = SHADER_PATH+"subface_vertex.glsl";
	  const std::string fragmentShader = SHADER_PATH+"subface_fragment.glsl";
	  INFO_LOG("vertex shader will be loaded from: "<<vertexShader);
	  INFO_LOG("surface shader will be loaded from: "<<fragmentShader);

	  CHECK_DIR_EXIST(SHADER_PATH);
	  CHECK_FILE_EXIST(vertexShader);
	  CHECK_FILE_EXIST(fragmentShader);

	  pSubFaceGLSL = pGLSLProgramObject(new GLSLProgramObject);
	  pSubFaceGLSL->attachVertexShader(vertexShader);
	  pSubFaceGLSL->attachFragmentShader(fragmentShader);
	  pSubFaceGLSL->link();
	}
  }
}
