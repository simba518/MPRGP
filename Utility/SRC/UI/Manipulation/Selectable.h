#ifndef _SELECTABLE_H_
#define _SELECTABLE_H_

#include <boost/shared_ptr.hpp>
#include <vector>
using std::vector;

namespace QGLVEXT{
  
  /**
   * @class Selectable the base class for selection.
   * 
   */
  class Selectable {
	
  public:
	// return the total element for selection, called in BaseViewer
	virtual int totalEleNum ()const = 0;

	// draw all of the elements for selection using opengl's drawWithNames,
	// called in BaseViewer.
	virtual void drawWithNames ()const = 0;

	// prepare datas before perform selection, called in BaseSelCtrl.
	virtual void prepareSelection(){
	  
	}

  };

  typedef boost::shared_ptr<Selectable > pSelectable; 
  typedef boost::shared_ptr<const Selectable > pSelectable_const; 
  
}//end of namespace

#endif /* _SELECTABLE_H_ */
