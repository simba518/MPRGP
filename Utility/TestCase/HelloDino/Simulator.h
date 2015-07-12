#ifndef _SIMULATOR_H_
#define _SIMULATOR_H_

#include <boost/shared_ptr.hpp>
#include <set>
#include <eigen3/Eigen/Dense>
#include <TetMesh.h>
using namespace std;
using namespace Eigen;
using namespace UTILITY;

namespace SIMULATOR{

  /**
   * @class Simulator interface for the solid simulation.
   * 
   */
  class Simulator{
	
  public:
	virtual string name()const{
	  return "un-known-simulator";
	}
	virtual bool init(const string filename) = 0;
	virtual void setVolMesh(pTetMesh_const tetMesh) = 0;
	virtual bool precompute(){return true;}
	virtual void reset(){
	  clearExtForce();
	  removeAllConNodes();
	}

	virtual void setConNodes(const vector<int> &con_nodes) = 0;
	virtual void setUc(const VectorXd &uc) = 0;
	virtual void removeAllConNodes() = 0;

	virtual void setExtForceOfNode(const int node_id,const double f[3]){
	  static VectorXd full_ext;
	  full_ext.resize(getFullDisp().size());
	  full_ext.setZero();
	  full_ext[node_id*3+0] = f[0];
	  full_ext[node_id*3+1] = f[1];
	  full_ext[node_id*3+2] = f[2];
	  this->setExtForce(full_ext);
	}
	virtual void setExtForce(const VectorXd &f_ext) = 0;
	virtual void clearExtForce(){
	  VectorXd full_ext(getFullDisp().size());
	  full_ext.setZero();
	  setExtForce(full_ext);
	}

	virtual bool forward() = 0;

	virtual const VectorXd &getFullDisp()const = 0;
	virtual const VectorXd &getFullVelocity()const = 0;
	virtual VectorXd &getFullDisp() = 0;
	virtual VectorXd &getFullVelocity() = 0;
	virtual bool computeElasticForce(const VectorXd &u,VectorXd &f)const{
	  return false;
	}

  };
  typedef boost::shared_ptr<Simulator> pSimulator;
  
}//end of namespace

#endif /*_SIMULATOR_H_*/
