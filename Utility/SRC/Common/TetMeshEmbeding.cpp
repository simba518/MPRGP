#include "AuxTools.h"
#include "TetMeshEmbeding.h"
using namespace UTILITY;

TetMeshEmbeding::TetMeshEmbeding(){
  _objMesh = pObjmesh(new Objmesh());
  _tetMesh = pTetMesh(new TetMesh());
}
TetMeshEmbeding::TetMeshEmbeding(pTetMesh tetMesh,pObjmesh objMesh){
  setTetMesh(tetMesh);
  setObjmesh(objMesh);
}
void TetMeshEmbeding::setTetMesh(pTetMesh tetMesh){
  assert(tetMesh);
  _tetMesh = tetMesh;
}
void TetMeshEmbeding::setObjmesh(pObjmesh objMesh){
  assert(objMesh);
  _objMesh = objMesh;
  _objRestVerts = _objMesh->getVerts();
}

pTetMesh TetMeshEmbeding::getTetMesh(){
  return _tetMesh;
}
pObjmesh TetMeshEmbeding::getObjMesh(){
  return _objMesh;
}
pTetMesh_const TetMeshEmbeding::getTetMesh()const{
  return _tetMesh;
}
pObjmesh_const TetMeshEmbeding::getObjMesh()const{
  return _objMesh;
}

void TetMeshEmbeding::buildInterpWeights(){
  assert(_tetMesh);
  assert(_objMesh);
  _tetMesh->buildInterpWeights(_objMesh->getVerts(),_nodes,_weights);
}
void TetMeshEmbeding::interpolate(const VectorXd& u){
  assert(_tetMesh);
  assert(_objMesh);
  _tetMesh->interpolate(_nodes,_weights,u,_uTarget);
  assert_eq(_objRestVerts.size(),_uTarget.size());
  _uTarget += _objRestVerts;
  _objMesh->setVerts(_uTarget);
}

bool TetMeshEmbeding::loadObjMesh(const string filename){
  assert(_objMesh);
  const bool succ = _objMesh->load(filename);
  _objRestVerts = _objMesh->getVerts();
  return succ;
}
bool TetMeshEmbeding::loadTetMesh(const string filename){
  assert(_tetMesh);
  return _tetMesh->load(filename);
}
	
bool TetMeshEmbeding::loadWeights(const string fname){
  
  INFILE(fin,fname);
  bool succ = fin.is_open();
	  
  const int numElementVertices = 4;
  const int numTargetLocations = _objMesh->getVertsNum();
  _nodes.resize(numElementVertices*numTargetLocations);
  _weights.resize(_nodes.size());

  int numVertices = -1;
  int currentVertex;

  // read the elements one by one and accumulate entries
  while (numVertices < numTargetLocations-1 && succ){
	numVertices++;
	if (fin.eof()){
	  ERROR_LOG("interpolation file("<<fname<<") is too short:"<<numVertices-1<<"<"<<numTargetLocations);
	  succ = false;
	}
	fin >> currentVertex;
	if (currentVertex != numVertices){
	  ERROR_LOG("consecutive vertex index at position: "<<currentVertex<<" mismatch.");
	  succ = false;
	}
	for(int j=0; j<numElementVertices; j++){
	  fin >> _nodes[currentVertex * numElementVertices+j];
	  fin >> _weights[currentVertex*numElementVertices+j];
	  ERROR_LOG_COND("failed to read data in file: "<<fname,fin.good());
	}
	succ &= fin.good();
  }

  fin.close();
  return succ;
}
bool TetMeshEmbeding::writeWeights(const string fname)const{

  OUTFILE(outf,fname);
  bool succ = outf.is_open();
  assert_eq((int)_nodes.size(),_weights.size());
  assert_eq(_nodes.size()%4,0);
  for (size_t i = 0; i < _nodes.size()/4 && succ; ++i){
	outf << i << " ";
	for (int j = 0; j < 4; ++j){
	  const int c = i*4+j;
	  outf << _nodes[c] << " " << _weights[c] <<" ";
	}
	outf << endl;
	succ = outf.good();
  }
  outf.close();
  return succ;
}

void TetMeshEmbeding::getBBox(double min[3],double max[3])const{
  assert(_tetMesh);
  assert(_objMesh);
  BBoxD b1 = _tetMesh->getBBox();
  const BBoxD b2 = _objMesh->getBBox();
  b1.add(b2);
  b1.getMaxConner(max);
  b1.getMinConner(min);
}
double TetMeshEmbeding::getMaxRadius()const{

  double min[3],max[3];
  getBBox(min,max);
  double r = max[0] - min[0];
  r = (r>=(max[1]-min[1]) ? r:(max[1]-min[1]));
  r = (r>=(max[2]-min[2]) ? r:(max[2]-min[2]));
  return r >= 0 ? r : 0;
}

const vector<int> &TetMeshEmbeding::getInterpNodes()const{
  return _nodes;
}
const VectorXd &TetMeshEmbeding::getInterpWeights()const{
  return _weights;
}
