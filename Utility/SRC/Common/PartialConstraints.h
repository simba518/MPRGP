#ifndef _PARTIALCONSTRAINTS_H_
#define _PARTIALCONSTRAINTS_H_

#include <boost/shared_ptr.hpp>
#include <eigen3/Eigen/Dense>
#include <SelectionGroup.h>
#include <algorithm> 
#include <fstream>
#include <boost/foreach.hpp>
#include <assertext.h>
using namespace std;
using namespace Eigen;

namespace UTILITY{
  
  /**
   * @class PartialConstraints partial constraints of one frame. Though the
   * constrained nodes are grouped, there are target positions for each node,
   * e.g not the barycenter constraints.
   * 
   */
  class PartialConstraints{
	
  public:
	PartialConstraints(const int frameid=0){
	  _frameid = frameid;
	}
	PartialConstraints(const vector<set<int> > &con_groups,const Matrix<double,3,-1> &pc,const int frame_id=0){
	  setConNodes(con_groups,pc,frame_id);
	}

	void setConNodes(const vector<set<int> > &con_groups,const Matrix<double,3,-1> &pc,const int frame_id=0){
	  _con_groups.setGroup( con_groups );
	  _pc = pc;
	  _frameid = frame_id;
	  assert_eq(_con_groups.numNodes(),_pc.cols());
	}
	void addConNodes(const vector<int> &con_group){
	  _con_groups.addGroup(con_group);
	}
	void rmConNodes(const vector<int> &con_group){
	  _con_groups.removeGroup(con_group);
	}
	void updatePc(const Matrix<double,3,-1> &pc){
	  assert_eq(pc.cols(),_con_groups.numNodes());
	  _pc = pc;
	}
	void updatePc(const Matrix<double,3,-1> &pc,const int group_id){

	  assert_in(group_id,0,numGroup()-1);
	  const int nodes = getConNodesSet()[group_id].size();
	  assert_eq(pc.cols(),nodes);
	  int c0 = 0;
	  for (int i = 0; i < group_id; ++i)
		c0 += getConNodesSet()[i].size();
	  _pc.block(0,c0,3,nodes) = pc;
	}

	const vector<set<int> > &getConNodesSet()const{
	  return _con_groups.getGroup();
	}
	const Matrix<double,3,-1> &getPc()const{
	  return _pc;
	}
	const Matrix<double,3,-1> getPc(const int group_id)const{
	  assert_in(group_id,0,numGroup()-1); 
	  const int nodes = getConNodesSet()[group_id].size();
	  int c0 = 0;
	  for (int i = 0; i < group_id; ++i)
		c0 += getConNodesSet()[i].size();
	  assert_in(c0,0,_pc.cols()-1);
	  assert_in(c0+nodes,0,_pc.cols());
	  return _pc.block(0,c0,3,nodes);
	}
	template<class VECTOR_INT, class VECTOR_DOUBLE>
	void getPartialCon(VECTOR_INT &conNodes, VECTOR_DOUBLE &conPos){
	  const int conNodesNum = _pc.cols();
	  if (conNodesNum > 0){
		_con_groups.getVector(conNodes);
		conPos.resize(conNodesNum*3);
		for (int c = 0; c < _pc.cols(); ++c){
		  conPos[c*3+0] = _pc.col(c)[0];
		  conPos[c*3+1] = _pc.col(c)[1];
		  conPos[c*3+2] = _pc.col(c)[2];
		}
	  }else{
		conNodes.clear();
		conPos.resize(0);
	  }
	}

	int numConNodes()const{
	  return _con_groups.numNodes();
	}
	int numGroup()const{
	  return getConNodesSet().size();
	}
	bool isEmpty()const{
	  return _con_groups.isEmpty();
	}
	const int getFrameId()const{
	  return _frameid;
	}

	void clear(){
	  _con_groups.clear();
	  _pc.resize(3,0);
	  _frameid = 0;
	}
	void setZero(){
	  _pc.setZero();
	}

	bool operator< (const PartialConstraints &other)const{
	  return this->getFrameId() < other.getFrameId();
	}

	bool write(const string filename)const;
	bool write(ofstream &outf)const;
	bool load(const string filename);
	bool load(ifstream &outf);

  private:
	int _frameid;
	SelectionGroup _con_groups;
	Matrix<double,3,-1> _pc;
  };
  typedef boost::shared_ptr<PartialConstraints> pPartialConstraints;
  typedef boost::shared_ptr<const PartialConstraints> pPartialConstraints_const;

  /**
   * @class ConNodesOfFrameSet a group of PartialConstraints, used to represent
   * the constraints of all the frames in animation editing.
   * 
   */
  class PartialConstraintsSet{
	
  public:
	void addPartialCon(const PartialConstraintsSet &parcons){
	  const set<pPartialConstraints> &cons = parcons.getPartialConSet();
	  BOOST_FOREACH(pPartialConstraints ele, cons){
		addPartialCon(ele);
	  }
	}
	void addPartialCon(const pPartialConstraints node_group){
	  removePartialCon(node_group->getFrameId());
	  _partialCons.insert(node_group);
	}
	void addConNodes(const vector<int> &con_nodes,const int frame_id){
	  pPartialConstraints pc = getPartialCon(frame_id);
	  if (!pc){
		pc = pPartialConstraints(new PartialConstraints(frame_id) );
		addPartialCon(pc);
	  }
	  pc->addConNodes(con_nodes);
	}
	void rmConNodes(const vector<int> &con_nodes,const int frame_id){
	  pPartialConstraints pc = getPartialCon(frame_id);
	  if (pc){
		pc->rmConNodes(con_nodes);
		if(pc->isEmpty())
		  removePartialCon(frame_id);
	  }
	}
	bool updatePc(const Matrix<double,3,-1> &pc,const int frame_id){
	  pPartialConstraints pcon = getPartialCon(frame_id);
	  if(pcon){
		pcon->updatePc(pc);
		return true;
	  }
	  return false;
	}

	pPartialConstraints_const getPartialCon(const int frame_id)const{
	  pPartialConstraints_const pc;
	  BOOST_FOREACH(pPartialConstraints_const ele, _partialCons)
		if ( ele->getFrameId() == frame_id){
		  pc = ele;
		  break;
		}
	  return pc;
	}
	pPartialConstraints getPartialCon(const int frame_id){
	  pPartialConstraints pc;
	  BOOST_FOREACH(pPartialConstraints ele, _partialCons)
		if ( ele->getFrameId() == frame_id){
		  pc = ele;
		  break;
		}
	  return pc;
	}
	const set<pPartialConstraints> &getPartialConSet()const{
	  return _partialCons;
	}
	set<pPartialConstraints> &getPartialConSet(){
	  return _partialCons;
	}
	vector<PartialConstraints> getPartialConSetSorted()const{
	  vector<PartialConstraints> conVec;
	  BOOST_FOREACH(pPartialConstraints con, _partialCons)
		if (!con->isEmpty())
		  conVec.push_back(*con);
	  std::sort (conVec.begin(), conVec.end());
	  return conVec;
	}
	void removePartialCon(const int frame_id){
	  BOOST_FOREACH(pPartialConstraints ele, _partialCons)
		if ( ele->getFrameId() == frame_id){
		  _partialCons.erase(ele);
		  break;
		}
	}

	void clear(){
	  _partialCons.clear();
	}
	
	bool load(const string filename);
	bool write(const string filename)const;

  private:
	set<pPartialConstraints> _partialCons;
  };
  
}//end of namespace

#endif /*_PARTIALCONSTRAINTS_H_*/
