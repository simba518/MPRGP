#define NV_REPORT_COMPILE_ERRORS

#include "GLSLProgramObject.h"
#include "NVShaderUtils.h"
using namespace UTILITY;

// ======================================================================================
GLSLProgramObject::GLSLProgramObject() :
  _progId(0)
{
}

// ======================================================================================
GLSLProgramObject::~GLSLProgramObject()
{
  destroy();
}

// ======================================================================================
void GLSLProgramObject::destroy()
{
  for (unsigned i = 0; i < _vertexShaders.size(); i++) {
	glDeleteShader(_vertexShaders[i]);
  }
  _vertexShaders.clear();
  _vertexShaders.resize( 0 );
  for (unsigned i = 0; i < _fragmentShaders.size(); i++) {
	glDeleteShader(_fragmentShaders[i]);
  }
  _fragmentShaders.clear();
  _fragmentShaders.resize( 0 );
  if (_progId != 0) {
	glDeleteProgram(_progId);
  }
}

// ======================================================================================
void GLSLProgramObject::attachVertexShader(std::string filename){

  std::string resolved_path;
  resolved_path = filename;
  GLuint shaderId = nv::CompileGLSLShaderFromFile(GL_VERTEX_SHADER, resolved_path.c_str());
  if (shaderId == 0) {
	std::cerr << "Error: Vertex shader failed to compile" << std::endl;
	exit(1);
  }
  _vertexShaders.push_back(shaderId);
}

// ======================================================================================
void GLSLProgramObject::attachFragmentShader(std::string filename){

  std::string resolved_path = filename;
  GLuint shaderId = nv::CompileGLSLShaderFromFile(GL_FRAGMENT_SHADER, resolved_path.c_str());
  if (shaderId == 0) {
	std::cerr << "Error: Fragment shader failed to compile" << std::endl;
	exit(1);
  }
  _fragmentShaders.push_back(shaderId);
}

// ======================================================================================
void GLSLProgramObject::link()
{
  _progId = glCreateProgram();

  for (unsigned i = 0; i < _vertexShaders.size(); i++) {
	glAttachShader(_progId, _vertexShaders[i]);
  }

  for (unsigned i = 0; i < _fragmentShaders.size(); i++) {
	glAttachShader(_progId, _fragmentShaders[i]);
  }

  glLinkProgram(_progId);

  GLint success = 0;
  glGetProgramiv(_progId, GL_LINK_STATUS, &success);

  if (!success) {
	char temp[1024];
	glGetProgramInfoLog(_progId, 1024, NULL, temp);
	printf("Failed to link program:\n%s\n", temp);
	exit(1);
  }
}

// ======================================================================================
void GLSLProgramObject::bind()
{
  glUseProgram(_progId);
}

// ======================================================================================
void GLSLProgramObject::unbind()
{
  glUseProgram(0);
}

// ======================================================================================
void GLSLProgramObject::setUniform(std::string name, GLfloat* val, int count)
{
  GLint id = glGetUniformLocation(_progId, name.c_str());
  if (id == -1) {
#ifdef NV_REPORT_UNIFORM_ERRORS
	std::cerr << m_vName << std::endl << m_fName << ":" << std::endl;
	std::cerr << "Warning: Invalid uniform parameter " << name << std::endl;
#endif
	return;
  }
  switch (count) {
  case 1:
	glUniform1fv(id, 1, val);
	break;
  case 2:
	glUniform2fv(id, 1, val);
	break;
  case 3:
	glUniform3fv(id, 1, val);
	break;
  case 4:
	glUniform4fv(id, 1, val);
	break;
  }
}

// ======================================================================================
void GLSLProgramObject::setTextureUnit(std::string texname, int texunit)
{
  GLint linked;
  glGetProgramiv(_progId, GL_LINK_STATUS, &linked);
  if (linked != GL_TRUE) {
	std::cerr << "Error: setTextureUnit needs program to be linked." << std::endl;
	exit(1);
  }
  GLint id = glGetUniformLocation(_progId, texname.c_str());
  if (id == -1) {
#ifdef NV_REPORT_UNIFORM_ERRORS
	std::cerr << "Warning: Invalid texture " << texname << std::endl;
#endif
	return;
  }
  glUniform1i(id, texunit);
}

// ======================================================================================
void GLSLProgramObject::bindTexture(GLenum target, std::string texname, GLuint texid, int texunit)
{
  glActiveTexture(GL_TEXTURE0 + texunit);
  glBindTexture(target, texid);
  setTextureUnit(texname, texunit);
  glActiveTexture(GL_TEXTURE0);
}
