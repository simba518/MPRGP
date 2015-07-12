#ifndef _SELFRENDERELE_H_
#define _SELFRENDERELE_H_

#ifdef WIN32
#define _interlockedbittestandset fk_m$_set
#define _interlockedbittestandreset fk_m$_reset
#include <boost/shared_ptr.hpp>
#undef _interlockedbittestandset
#undef _interlockedbittestandreset
#else
#include <boost/shared_ptr.hpp>
#endif

namespace QGLVEXT{

  /**
   * The render priority of the object. Object with priority of FIRST_RENDER
   * will be rendered first, and those with FOURTH_RENDER will be rendered last.
   * Objects with similar priority will be rendered similar to the oders of
   * being regested to the viewer. The default rendered priority is
   * SECOND_RENDER.
   * 
   * @note the render order is compressponding to QGLViewerExt::selfRender(), so
   * if new order is added, implementation of that function should be changed.
   * 
   */
  enum RENDER_PRIORITY {FIRST_RENDER, SECOND_RENDER, THIRD_RENDER, FOURTH_RENDER};

  /**
   * @class SelfRenderEle the base class for rendering. all the rendering
   * code should be warped in the child class. Then for rendering, the child
   * object should be registered to the viewer (using
   * BaseViewer::addSelfRenderEle(..)).
   * 
   * This class also inherent from the class BaseTextRecorder, and providing
   * interfaces for recording text for display.
   * @see BaseViewer
   * @see BaseTextRecorder
   */
  class SelfRenderEle{
	
  public:
	SelfRenderEle(const RENDER_PRIORITY render_priority = SECOND_RENDER){
	  setRenderPriority(render_priority);
	}

	void setRenderPriority(const RENDER_PRIORITY render_priority){
	  this->render_priority = render_priority;
	}

	RENDER_PRIORITY getRenderPriority()const{
	  return render_priority;
	}

	virtual void draw()const = 0;

	// non-const version of draw
	virtual void render(){
	  this->draw();
	}
	
  private:
	RENDER_PRIORITY render_priority;

  };

  typedef boost::shared_ptr< const SelfRenderEle > pSelfRenderEle_const; 
  typedef boost::shared_ptr< SelfRenderEle > pSelfRenderEle; 

}//end of namespace

#endif /* _SELFRENDERELE_H_ */
