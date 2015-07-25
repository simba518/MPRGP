#ifndef _TETMESHEMBEDING_H_
#define _TETMESHEMBEDING_H_

#include <TetMesh.h>
#include <Objmesh.h>

namespace UTILITY{

  class TetMeshEmbeding{

  public:
	TetMeshEmbeding();
	TetMeshEmbeding(pTetMesh tetMesh,pObjmesh objMesh);
	void setTetMesh(pTetMesh tetMesh);
	void setObjmesh(pObjmesh objMesh);

	pTetMesh getTetMesh();
	pObjmesh getObjMesh();
	pTetMesh_const getTetMesh()const;
	pObjmesh_const getObjMesh()const;

	void buildInterpWeights();
	void interpolate(const VectorXd& u);

	bool loadObjMesh(const string filename);
	bool loadTetMesh(const string filename);
	
	bool loadWeights(const string fname);
	bool writeWeights(const string fname)const;

	void getBBox(double min[3],double max[3])const;
	double getMaxRadius()const;

	const vector<int> &getInterpNodes()const;
	const VectorXd &getInterpWeights()const;

  private:
	pTetMesh _tetMesh;
	pObjmesh _objMesh;
	VectorXd _objRestVerts;

	vector<int> _nodes;
	VectorXd _weights;

	VectorXd _uTarget;
  };
  
  typedef boost::shared_ptr<TetMeshEmbeding> pTetMeshEmbeding;
  typedef boost::shared_ptr<const TetMeshEmbeding> pTetMeshEmbeding_const;

}//end of namespace

#endif /* _TETMESHEMBEDING_H_ */
