#ifndef _GLSLPROGRAMOBJECT_H_
#define _GLSLPROGRAMOBJECT_H_

#include <boost/shared_ptr.hpp>
#include <GL/glew.h>
#include <string>
#include <iostream>
#include <vector>

namespace UTILITY{

  class GLSLProgramObject{

  public:
	GLSLProgramObject();
	virtual ~GLSLProgramObject();

	void destroy();
	void bind();
	void unbind();

	void setUniform(std::string name, GLfloat* val,int count);
	void setTextureUnit(std::string texname, int texunit);
	void bindTexture(GLenum target,std::string texname,GLuint texid,int texunit);

	void bindTexture2D(std::string texname, GLuint texid, int texunit) {
	  bindTexture(GL_TEXTURE_2D, texname, texid, texunit);
	}

	void bindTexture3D(std::string texname, GLuint texid, int texunit) {
	  bindTexture(GL_TEXTURE_3D, texname, texid, texunit);
	}

	void bindTextureRECT(std::string texname, GLuint texid, int texunit) {
	  bindTexture(GL_TEXTURE_RECTANGLE_ARB, texname, texid, texunit);
	}

	void attachVertexShader(std::string filename);

	void attachFragmentShader(std::string filename);

	void link();

	inline GLuint getProgId() { return _progId; }
	
  protected:
	std::vector<GLuint>		_vertexShaders;
	std::vector<GLuint>		_fragmentShaders;
	GLuint					_progId;
  };

  typedef boost::shared_ptr<GLSLProgramObject> pGLSLProgramObject;
}

#endif /* _GLSLPROGRAMOBJECT_H_ */
