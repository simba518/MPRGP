#ifndef _PHONGSHADER_H_
#define _PHONGSHADER_H_

#include <string>
#include <boost/shared_ptr.hpp>
#include <GLSLProgramObject.h>

namespace UTILITY{
  
  class PhongShader{
	
  public:
	static boost::shared_ptr<PhongShader> getInstance(){
	  
	  if(p_instance == NULL){
		p_instance = boost::shared_ptr<PhongShader>(new PhongShader());
	  }
	  return p_instance;
	}
	/// @note this function should be called after OpenGL initialized.
	void initialize();
	void bind(){
	  if (pSubFaceGLSL)
		pSubFaceGLSL->bind();
	}
	void unbind(){
	  if (pSubFaceGLSL)
		pSubFaceGLSL->unbind();
	}
	
  protected:
    PhongShader(){}
	
  private:
	pGLSLProgramObject pSubFaceGLSL;
	static boost::shared_ptr<PhongShader> p_instance;
  };
  
  typedef boost::shared_ptr<PhongShader> pPhongShader;
  
}//end of namespace

#endif /*_PHONGSHADER_H_*/
