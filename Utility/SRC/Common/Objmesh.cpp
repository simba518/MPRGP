#include <boost/filesystem.hpp>
#include <Objmesh.h>
#include <AuxTools.h>
using boost::filesystem::path;
using namespace UTILITY;
using namespace Eigen;

void ObjMtl::setDefault(){

  name = "defualt";

  ambient[0] = 0.0f;       ambient[1] = 0.1f;       ambient[2] = 0.0f;
  diffuse[0] = 0.0f;       diffuse[1] = 0.6f;       diffuse[2] = 0.0f;
  specular[0] = 0.35f;      specular[1] = 0.35f;      specular[2] = 0.35f;
  transmittance[0] = 1.0f; transmittance[1] = 1.0f; transmittance[2] = 1.0f;
  emission[0] = 0.0f;      emission[1] = 0.0f;      emission[2] = 0.0f;

  shininess = 200;
  ior = 1;
}

// v : 0
// vn: 1
// vt: 2
// f : 3
// other: -1
int Objmesh::type(const string &line)const{
  
  int t = -1;
  if(line.size() > 2){
	if('v' == line[0]){
	  t = 0;
	  if('n' == line[1])
		t = 1;
	  else if('t' == line[1])
		t = 2;
	}else if('f' == line[0])
	  t = 3;
  }
  return t;
}

// return value
// 1: 1 2 3
// 2: 1//1 2//2 3//3
// 3: 1/1/1 2/1/2 3/3/3
int getFaceNormal(string &line,vector<int>&faceNormal){

  int faceType = 1;
  faceNormal.clear();
  for (size_t i = 0; i < line.size(); ++i){
    if('/' == line[i]){
	  faceType = 3;
	  line[i] = ' ';
	  if('/' == line[i+1]){
		faceType = 2;
		line[i+1] = ' ';
	  }
	}
  }
  istringstream ss(line);
  string temp;
  ss >> temp;
  int i;
  while(!ss.eof()&&ss.good()){
	ss >> i;
	faceNormal.push_back(i-1);
  }
  return faceType;
}

// load triangle obj file
// 1. load all meshes as triangle mesh.
// 2. load only one mtl.
// 3. no textures.
bool Objmesh::load(const string fname){

  std::ifstream inf;
  inf.open(string(fname).c_str());
  if(!inf.is_open()){
	ERROR_LOG("failed to open file for output: "<<fname); 
	return false;  
  }
  
  // count the number of v,vn,vt and f.
  int numType[4] = {0}; // num of v,vn,vt,f
  string line;
  string mtlfile;
  while(!inf.eof()&&getline(inf,line)){
	const int t = type(line);
	if(t >=0 && t < 4){
	  numType[t]++;
	}else if(line.find("mtllib") != std::string::npos){
	  const string temp = path(fname).parent_path().string();
	  const string dir = temp+"//";
	  mtlfile = dir + line.substr(7,line.size()-7);
	}
  }

  vector<double> verts(numType[0]*3,0);
  vector<double> vertNormal(numType[1]*3,0);
  vector<int> faces, normalIndex;
  faces.reserve(numType[3]*3);
  normalIndex.reserve(numType[3]*3);

  // load real data
  inf.clear();
  inf.seekg(0,ios::beg);
  int v=0,vn=0;
  bool succ = true;
  while(!inf.eof()&&getline(inf,line) && succ){

	const int t = type(line);
	if(t >= 0){
	  istringstream ss(line);
	  string temp; ss >> temp;
	  switch (t){
	  case 0:{
		ss >> verts[v] >> verts[v+1] >> verts[v+2];
		v += 3;
		break;
	  };
	  case 1:{
		ss >> vertNormal[vn] >> vertNormal[vn+1] >> vertNormal[vn+2];
		vn += 3;
		break;
	  };
	  case 3:{
		static vector<int> faceNormal;
		const int ft = getFaceNormal(line,faceNormal);
		assert_gt(ft,0);
		for (size_t i = 1; i < faceNormal.size()/ft-1; i++){
		  assert_in(i*ft+ft,0,faceNormal.size());
		  const int p = i*ft;
		  const int q = p+ft;
		  faces.push_back(faceNormal[0]);
		  faces.push_back(faceNormal[p]);
		  faces.push_back(faceNormal[q]);
		  if (ft >= 2){
			assert_in(q+ft-1,0,faceNormal.size());
			normalIndex.push_back(faceNormal[ft-1]);
			normalIndex.push_back(faceNormal[p+ft-1]);
			normalIndex.push_back(faceNormal[q+ft-1]);
		  }
		}
		break;
	  };
	  }
	}
  }
  inf.close();

  setVerts(verts);
  setVertNormals(vertNormal);
  setFaces(faces);
  setNormalIndex(normalIndex);

  // load mtl
  if(mtlfile.size() > 0){
	if(!loadMtl(mtlfile)){
	  _mtl.setDefault();
	  WARN_LOG("failed to load the mtl from "<<mtlfile<<"default mtl will be used.");
	}
  }

  return succ;
}

bool Objmesh::write(const string fname)const{

  ofstream outf;
  outf.open(fname.c_str());
  if(!outf.is_open()){
	ERROR_LOG("failed to open file "<< fname <<" for writing.");
	return false;  
  }

  const int nv = getVertsNum();
  const int nf = getFacesNum();
  outf << "# number of vertices "<< nv <<endl;
  outf << "# number of faces "<< nf<< endl;
  if(_vertNormal.size() == nv*3)
	outf << "# number of normals " << nv << endl;
  outf << endl;

  for (int i = 0; i < nv; ++i){
	outf<< "v " << _verts[i*3] << " ";
	outf<< _verts[i*3+1] << " ";
	outf<< _verts[i*3+2] << endl;
  }
  outf << endl;

  if(_vertNormal.size() == nv*3){
	for (int i = 0; i < nv; ++i){
	  outf<< "vn " << _vertNormal[i*3] << " ";
	  outf<< _vertNormal[i*3+1] << " ";
	  outf<< _vertNormal[i*3+2] << endl;
	}
	outf << endl;
	for (int i = 0; i < nf; ++i){
	  outf<<"f "<< _faces[i*3]+1<<"//"<<_faces[i*3]+1<<" ";
	  outf<< _faces[i*3+1]+1<<"//"<<_faces[i*3+1]+1<<" ";
	  outf<< _faces[i*3+2]+1<<"//"<<_faces[i*3+2]+1<<endl;
	}
  }else{
	for (int i = 0; i < nf; ++i){
	  outf<<"f "<< _faces[i*3]+1<<" ";
	  outf<< _faces[i*3+1]+1<<" ";
	  outf<< _faces[i*3+2]+1<<endl;
	}
  }

  outf.close();
  return true;
}

bool Objmesh::writeVTK(const std::string& filename,const bool binary)const{

  VTKWriter writer(binary);
  writer.addPoints(_verts);
  writer.addTriangles(_faces);
  return writer.write(filename);
}

bool Objmesh::loadMtl(const string fname){
  
  std::ifstream inf;
  inf.open(string(fname).c_str());
  if(!inf.is_open()){
	ERROR_LOG("failed to open file for output: "<<fname);
	return false;  
  }
  string key;
  while(!inf.eof()&&(inf >> key).good()){
	if(key.size()<=0||'#' == key[0]){
	  continue;
	}else if(string("Kd") == key){
	  inf >> _mtl.diffuse[0] >> _mtl.diffuse[1] >> _mtl.diffuse[2];
	}else if(string("Ka") == key){
	  inf >> _mtl.ambient[0] >> _mtl.ambient[1] >> _mtl.ambient[2];
	}else if(string("Tf") == key){
	  inf >> _mtl.transmittance[0] >> _mtl.transmittance[1] >> _mtl.transmittance[2];
	}else if(string("Ks") == key){
	  inf >> _mtl.specular[0] >> _mtl.specular[1] >> _mtl.specular[2];
	}else if(string("Ke") == key){
	  inf >> _mtl.emission[0] >> _mtl.emission[1] >> _mtl.emission[2];
	}else if(string("Ns") == key){
	  inf >> _mtl.shininess;
	}else if(string("Ni") == key){
	  inf >> _mtl.ior;
	}else if(string("newmtl") == key){
	  inf >> _mtl.name;
	}
  }
  inf.close();
  return true;
}

void Objmesh::moveCenterTo(const Eigen::Vector3d &pos){
  
  if (_verts.size() > 0){
	const Eigen::Vector3d disp = pos-getCenter();
	for (int i = 0; i < _verts.size()/3; ++i){
	  _verts.segment<3>(i*3) += disp;
	}
  }
}

void Objmesh::combine(const Objmesh &other){

  const VectorXd &v = other.getVerts();
  const VectorXd &n = other.getVertNormal();
  const VectorXi &fi = other.getFaces();
  const VectorXi &ni = other.getNormalIndex();
  
  VectorXd t_verts(v.size()+_verts.size()), t_vertNormal(n.size()+_vertNormal.size());
  VectorXi t_faces(fi.size()+_faces.size()), t_normalIndex(ni.size()+_normalIndex.size());

  t_verts.head(_verts.size()) = _verts;
  t_verts.tail(v.size()) = v;

  t_vertNormal.head(_vertNormal.size()) = _vertNormal;
  t_vertNormal.tail(n.size()) = n;

  t_faces.head(_faces.size()) = _faces;
  t_faces.tail(fi.size()).setConstant(_verts.size()/3);
  t_faces.tail(fi.size()) += fi;

  t_normalIndex.head(_normalIndex.size()) = _normalIndex;
  t_normalIndex.tail(ni.size()).setConstant(_verts.size()/3);
  t_normalIndex.tail(ni.size()) += ni;

  _verts = t_verts;
  _vertNormal = t_vertNormal;
  _faces = t_faces;
  _normalIndex = t_normalIndex;
}

void Objmesh::toPlanes(const Eigen::VectorXd &verts, const Eigen::VectorXi &faces, VVec4d &planes){
 
  assert_eq(verts.size()%3,0);
  assert_eq(faces.size()%3,0);

  planes.resize(faces.size()/3);

  for (int f = 0; f < faces.size()/3; ++f){

	const int i0 = faces[f*3+0];
	const int i1 = faces[f*3+1];
	const int i2 = faces[f*3+2];

	assert_in(i0*3,0,verts.size()-3);
	assert_in(i1*3,0,verts.size()-3);
	assert_in(i2*3,0,verts.size()-3);

	const Vector3d v0 = verts.segment<3>(i0*3);
	const Vector3d v1 = verts.segment<3>(i1*3);
	const Vector3d v2 = verts.segment<3>(i2*3);

	planes[f].head(3) = (v1-v0).cross(v2-v0).normalized();
	planes[f][3] = -v0.dot(planes[f].head(3));

  }
}
