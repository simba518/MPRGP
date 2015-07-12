#ifndef _OBJMESH_H_
#define _OBJMESH_H_

#include <map>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <Eigen/Dense>
#include <BBox.h>
#include <assertext.h>
#include <VTKWriter.h>
using namespace std;

namespace UTILITY{

  typedef std::vector<Eigen::Vector4d,Eigen::aligned_allocator<Eigen::Vector4d> > VVec4d;

  typedef struct _ObjMtl{

	_ObjMtl(){
	  setDefault();
	}
	void setDefault();

	string name;

    float ambient[3];
    float diffuse[3];
    float specular[3];
    float transmittance[3];
    float emission[3];
    float shininess;
    float ior;  // index of refraction

	string ambient_texname;
    string diffuse_texname;
    string specular_texname;
    string normal_texname;
	
  }ObjMtl;

  class Objmesh{
	
  public:
	Objmesh(){}
	template<class VECTOR_D, class VECTOR_I>
	Objmesh(const VECTOR_D&vertices,const VECTOR_I&faces){
	  assert_eq(vertices.size()%3,0);
	  setVerts(vertices);
	  setFaces(faces);
	}
	
	// set
	void setMtl(const ObjMtl &mtl){
	  _mtl = mtl;
	}
	void setVerts(const Eigen::VectorXd &v){ _verts = v;}
	template<class VECTOR>
	void setVerts(const VECTOR &v){setVec(v,_verts);}
	template<class VECTOR>
	void setVertNormals(const VECTOR&vn){setVec(vn,_vertNormal);}
	template<class VECTOR>
	void setFaces(const VECTOR &f){setVec(f,_faces);}
	template<class VECTOR>
	void setNormalIndex(const VECTOR &ni){setVec(ni,_normalIndex);}
	void combine(const Objmesh &other); // add the data from other to this.

	// get
	const ObjMtl &getMtl()const{
	  return _mtl;
	}
	const Eigen::VectorXd &getVerts()const{
	  return _verts;
	}
	const Eigen::VectorXd &getVertNormal()const{
	  return _vertNormal;
	}
	const Eigen::VectorXi &getFaces()const{
	  return _faces;
	}
	const Eigen::VectorXi &getNormalIndex()const{
	  return _normalIndex;
	}
	const Eigen::Vector3d getVerts(const int i)const{
	  return getSubV3(_verts,i);
	}
	const Eigen::Vector3d getVertNormal(const int i)const{
	  return getSubV3(_vertNormal,i);
	}
	const Eigen::Vector3i getFaces(const int i)const{
	  return getSubV3(_faces,i);
	}
	const Eigen::Vector3i getNormalIndex(const int i)const{
	  return getSubV3(_normalIndex,i);
	}
	int getVertsNum()const{
	  return _verts.size()/3;
	}
	int getFacesNum()const{
	  return _faces.size()/3;
	}
	int getNormalNum()const{
	  return _vertNormal.size()/3;
	}

	BBoxD getBBox()const{
	  BBoxD box;
	  if(_verts.size()>0){
		assert_eq(_verts.size()%3,0);
		box.reset(&(_verts[0]),_verts.size()/3);
	  }
	  return box;
	}
	Eigen::Vector3d getCenter()const{
	  Eigen::Vector3d center;
	  center.setZero();
	  for (int i = 0; i < getVertsNum(); ++i)  center += getVerts(i);
	  if(getVertsNum()>0)	center = center/getVertsNum();
	  return center;
	}
	void moveCenterTo(const Eigen::Vector3d &pos);
	void toPlanes(VVec4d &planes)const{
	  toPlanes(_verts,_faces,planes);
	}
	// convert to the planes: n[i]*x + p[i] = 0, where n[i]=planes[i][0,1,2], p[i] = planes[i][3], and n[i] is normalized.
	static void toPlanes(const Eigen::VectorXd &verts, const Eigen::VectorXi &faces, VVec4d &planes);

	// IO
	// load triangle obj file
	// 1. load mesh as triangle mesh.
	// 2. load only vertex normal.
	// 3. load only one mtl.
	// 4. no textures.
	bool load(const string fname);
	bool loadMtl(const string fname);
	bool write(const string fname)const;
	bool writeVTK(const std::string& filename,const bool binary=true)const;

  protected:
	template<class T> 
	Eigen::Matrix<T,3,1> getSubV3(const Eigen::Matrix<T,-1,1> &V,int i)const{
	  assert_in(i,0,V.size()/3-1);
	  return V.segment(i*3,3);
	}
	template<class VECTOR, class T> 
	void setVec(const VECTOR &inputV, Eigen::Matrix<T,-1,1> &V)const{
	  V.resize(inputV.size());
	  for (int i = 0; i < V.size(); ++i) V[i] = inputV[i];
	}
	int type(const string &line)const;
	
  private:
	Eigen::VectorXd _verts;
	Eigen::VectorXd _vertNormal;
	Eigen::VectorXi _faces;
	Eigen::VectorXi _normalIndex;
	ObjMtl _mtl;
  };
  
  typedef boost::shared_ptr<Objmesh> pObjmesh;
  typedef boost::shared_ptr<const Objmesh> pObjmesh_const;

}

#endif /*_OBJMESH_H_*/
