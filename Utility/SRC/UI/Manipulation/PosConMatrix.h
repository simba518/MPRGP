#ifndef _POSCONMATRIX_H_
#define _POSCONMATRIX_H_

#include <boost/foreach.hpp>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
using namespace Eigen;
#include <AssertExt.h>
#include <SelectionGroup.h>

namespace QGLVEXT{
  
  /**
   * @class PosConMatrix Given constraint nodes and current position of all the
   * nodes, assemble the sparse matrix C and target constraint positions xc.
   * 
   */
  template <typename T>
  class PosConMatrix{
	
  public:
	PosConMatrix(){}
	void setConstraints(const vector<set<int> > &groups, const Matrix<T,-1,1> &x){
	  _con_node_groups.setGroup(groups);
	  buildConMatrix(x.size());
	  buildConPositions(x);
	}
	void addConstraints(const vector<int> &nodes, const Matrix<T,-1,1> &x){
	  _con_node_groups.addGroup(nodes);
	  buildConMatrix(x.size());
	  buildConPositions(x);	
	}
	void rmConstraints(const vector<int> &nodes, const Matrix<T,-1,1> &x){
	  _con_node_groups.removeGroup(nodes);
	  buildConMatrix(x.size());
	  buildConPositions(x);
	}
	const SparseMatrix<int> &C()const{ return _C;}
	const Matrix<T,-1,1> &xc()const{ return _xc;}

  protected:
	void buildConMatrix(const int x_lengh){

	  if (_con_node_groups.numNodes()<=0){
		_C.resize(0,0);
		return ;
	  }

	  const vector<set<int> > &g = _con_node_groups.getGroup();
	  typedef Eigen::Triplet<int> Tri;
	  vector<Tri> triplets;
	  triplets.reserve(_con_node_groups.numNodes()*3);
	  int node_row = 0;
	  BOOST_FOREACH(const set<int> &s, g){
	  	BOOST_FOREACH(const int i, s){
	  	  ASSERTIN (i*3,0,x_lengh-3);
	  	  triplets.push_back(Tri(node_row*3+0, i*3+0, 1));
	  	  triplets.push_back(Tri(node_row*3+1, i*3+1, 1));
	  	  triplets.push_back(Tri(node_row*3+2, i*3+2, 1));
	  	  node_row ++;
	  	}
	  }
	  _C.resize(_con_node_groups.numNodes()*3, x_lengh);
	  _C.setZero();
	  _C.setFromTriplets(triplets.begin(), triplets.end());
	  _C.makeCompressed();
	}

	void buildConPositions(const Matrix<T,-1,1> &x){

	  if (_con_node_groups.numNodes()<=0){
		_xc.resize(0);
		return ;
	  }

	  const vector<set<int> > &g = _con_node_groups.getGroup();
	  int xci = 0;
	  _xc.resize(_con_node_groups.numNodes()*3);
	  BOOST_FOREACH(const set<int> &s, g){
		BOOST_FOREACH(const int i, s){
		  assert_in (i*3,0,x.size()-3);
		  _xc[xci*3+0] = x[i*3+0];
		  _xc[xci*3+1] = x[i*3+1];
		  _xc[xci*3+2] = x[i*3+2];
		  xci++;
		}
	  }
	}
	
  private:
	SelectionGroup _con_node_groups;
	SparseMatrix<int> _C;
	Matrix<T,-1,1> _xc;
  };
  
}//end of namespace

#endif /*_POSCONMATRIX_H_*/
