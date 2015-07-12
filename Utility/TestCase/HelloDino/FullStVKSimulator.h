#ifndef _FULLSTVKSIMULATOR_H_
#define _FULLSTVKSIMULATOR_H_

#include "Simulator.h"
#include <ConMatrixTools.h>
#include <FullSimulator.h>
using namespace UTILITY;

namespace SIMULATOR{
  
  /**
   * @class FullStVKSimulator solid simulation using full stvk model.
   * 
   */
  class FullStVKSimulator:public Simulator{
	
  public:
	FullStVKSimulator():sim_name("full stvk"){
	  stvkModel = pFullStVKSimModel(new FullStVKSimModel());
	  simulator = pBaseFullSim(new PenSemiImpFullSim(stvkModel));
	}
	FullStVKSimulator(pBaseFullSim sim,const string name):sim_name(name){
	  stvkModel = pFullStVKSimModel(new FullStVKSimModel());
	  simulator = sim;
	  simulator->setElasticModel(stvkModel);
	}
	string name()const{return sim_name;}
	bool init(const string filename){return simulator->init(filename);}
	void setVolMesh(pTetMesh_const tetMesh){stvkModel->setTetMesh(tetMesh);}
	bool precompute(){
	  simulator->prepare();
	  setConNodes(con_nodes);
	  VectorXd uc(con_nodes.size()*3);
	  uc.setZero();
	  setUc(uc);
	  return true;
	}
	void reset(){
	  clearExtForce();
	  removeAllConNodes();
	  simulator->reset();
	}

	void setConNodes(const vector<int> &nodes){

	  con_nodes = nodes;
	  if (nodes.size() <= 0){
		removeAllConNodes();
		return;
	  }

	  const int n = stvkModel->dimension()/3;
	  if (n > 0){
		VecT trip_C;
		UTILITY::computeConM(nodes, trip_C, n);
		simulator->setConM(trip_C, nodes.size()*3, n*3);
	  }
	}
	void setUc(const VectorXd &uc){simulator->setUc(uc);}
	void removeAllConNodes(){simulator->removeAllCon();}
	
	void setExtForceOfNode(const int nodeId,const double f[3]){
	  simulator->setExtForceOfNode(f, nodeId);
	}
	void setExtForce(const VectorXd &f_ext){simulator->setExtForce(f_ext);}
	void clearExtForce(){simulator->setExtForceForAllNodes(0.0f,0.0f,0.0f);}

	bool forward(){return simulator->forward();}

	const VectorXd &getFullDisp()const{return simulator->getU();}
	const VectorXd &getFullVelocity()const{
	  return simulator->getV();
	}
	VectorXd &getFullDisp(){return simulator->getU();}
	VectorXd &getFullVelocity(){
	  return simulator->getV();
	}
	bool computeElasticForce(const VectorXd &u,VectorXd &f)const{
	  if (stvkModel)
		return stvkModel->evaluateF(u,f);
	  return false;
	}
	
  private:
	vector<int> con_nodes;
	pFullStVKSimModel stvkModel;
	pBaseFullSim simulator;
	const string sim_name;
  };
  typedef boost::shared_ptr<FullStVKSimulator> pFullStVKSimulator;
  
}//end of namespace

#endif /*_FULLSTVKSIMULATOR_H_*/
