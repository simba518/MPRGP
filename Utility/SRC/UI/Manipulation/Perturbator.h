#ifndef _PERTURBATOR_H_
#define _PERTURBATOR_H_

#include <assert.h>
#include <stdlib.h>
#include <memory.h>

namespace QGLVEXT{

  /** 
   * @class Perturbator
   * @brief This class is used to produce the perturb force when user 
   * draging one of the vertex of volumetric mesh.
   * The draging forces is a 3x1 vertor which is computed as
   * F = (dragTo - start_pos)*drag_compilance.
   * where
   * both dragTo and start_pos are computed by unproject the mouse position to the world coord.
   * 1. dragTo is current position of the mouse.
   * 2. start_pos is the position when the perturb start.
   * 3. drag_compilance is a positive varialbe to scale the perturb forces.
   */
  class Perturbator{
	
  public:
	Perturbator(){
	  perturb_compilance = 1000.0f;
	  clear();
	}
	bool isDrag()const{
	  return (perturbNodeId >= 0);
	}
	void setPerturCompilance(double c){
	  assert(c > 0.0f);
	  perturb_compilance = c;
	}
	void setPerturbNodeId(int volmesh_node_id){
	  assert(volmesh_node_id >= 0);
	  perturbNodeId = volmesh_node_id;
	}
	void setStartPosition(double x,double y,double z){
	  start_pos[0] = x;
	  start_pos[1] = y;
	  start_pos[2] = z;
	  perturb_to[0] = x;
	  perturb_to[1] = y;
	  perturb_to[2] = z;
	}
	void dragTo(double x,double y,double z){
	  perturb_to[0] = x;
	  perturb_to[1] = y;
	  perturb_to[2] = z;	  
	}

	int getPerturbNodeId()const{
	  return perturbNodeId;
	}
	void getPerturbForce(double force[3])const{
	  assert(isDrag());
	  force[0] = perturb_compilance*(perturb_to[0] - start_pos[0]);
	  force[1] = perturb_compilance*(perturb_to[1] - start_pos[1]);
	  force[2] = perturb_compilance*(perturb_to[2] - start_pos[2]);
	}

	void stopDrag(){
	  clear();
	}
	void clear(){
	  perturbNodeId = -1;
	  memset((void*)start_pos,0,sizeof(double)*3);
	  memset((void*)perturb_to,0,sizeof(double)*3);
	}
	
  private:
	int perturbNodeId;
	double start_pos[3];
	double perturb_to[3];
	double perturb_compilance;
	
  };
  
}//end of namespace

#endif /*_PERTURBATOR_H_*/
