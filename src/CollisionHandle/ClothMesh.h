#ifndef CLOTH_MESH_H
#define CLOTH_MESH_H

#include <ObjMesh.h>
#include <Eigen/Sparse>
#include <eigen3/Eigen/Dense>
using namespace Eigen;

USE_PRJ_NAMESPACE

class ClothMesh : public Serializable
{
public:
    typedef Eigen::Matrix<scalarD,-1,1> Vec;
    enum MASS_MODE {
        BARYCENTER		=1,
        CIRCUMCENTER	=2,
    };
    enum MESH_TYPE {
        CLOTH_MESH=0,
        RIGID_MESH=1,
    };
    struct ClothTriangle;
    struct ClothVertex : public Serializable 
	{
		ClothVertex();
        ClothVertex(const Vec3d& pos,MESH_TYPE type);
		bool write(std::ostream& os,IOData* dat) const;
		bool read(std::istream& is,IOData* dat);
		boost::shared_ptr<Serializable> copy() const {
			return boost::shared_ptr<Serializable>(new ClothVertex);
		}
        //data
        Vec3d _pos,_lastPos,_vel,_pos0;
        scalarD _mass,_weight;
	    sizeType _index, body_index, node_id_on_body;
        //for collision detect
        sizeType _type;
        vector<boost::shared_ptr<ClothTriangle> > _oneRing;
    };
    struct ClothEdge : public Serializable 
	{
		ClothEdge();
        ClothEdge(boost::shared_ptr<ClothVertex> v0,boost::shared_ptr<ClothVertex> v1,MESH_TYPE type);
		bool write(std::ostream& os,IOData* dat) const;
		bool read(std::istream& is,IOData* dat);
		boost::shared_ptr<Serializable> copy() const {
			return boost::shared_ptr<Serializable>(new ClothEdge);
		}
	    Vec3d getNormal()const{
		  Vec3d n1 = _t[0]->getNormal();
		  Vec3d n2 = _t[1]->getNormal();
		  return (n1+n2)*0.5f;
		}
        //data
        boost::shared_ptr<ClothVertex> _v[2];
        boost::shared_ptr<ClothTriangle> _t[2];
        Vec3d _pos,_vel,_pos0;
        scalarD _mass;
        sizeType _index;
        //for collision detection
        sizeType _type;
    };
    struct ClothTriangle : public Serializable 
	{
		ClothTriangle();
        ClothTriangle(boost::shared_ptr<ClothEdge> e0,bool p0,boost::shared_ptr<ClothEdge> e1,bool p1,boost::shared_ptr<ClothEdge> e2,bool p2,MESH_TYPE type);
        boost::shared_ptr<ClothVertex> getV0() const;
        boost::shared_ptr<ClothVertex> getV1() const;
        boost::shared_ptr<ClothVertex> getV2() const;
		bool write(std::ostream& os,IOData* dat) const;
		bool read(std::istream& is,IOData* dat);
		boost::shared_ptr<Serializable> copy() const {
			return boost::shared_ptr<Serializable>(new ClothTriangle);
		}
		template <typename T>
		boost::shared_ptr<T> getElem(sizeType i) const{return _e[i];}
	    Vec3d getNormal()const{
		  const Vec3d &p0 = getV0()->_lastPos;
		  const Vec3d &p1 = getV1()->_lastPos;
		  const Vec3d &p2 = getV2()->_lastPos;
		  return (p1-p0).cross(p2-p0);
		}
		// template <>
		// boost::shared_ptr<ClothVertex> getElem(sizeType i) const{
		// 	if(i == 0)
		// 		return getV0();
		// 	else if(i == 1)
		// 		return getV1();
		// 	else return getV2();
		// }
        //data
        boost::shared_ptr<ClothEdge> _e[3];
        Vec3c _edgeDir;
        sizeType _index;
        //for collision detection
        sizeType _type;
		vector<bool>* _activeTag;
    };
public:
    ClothMesh();
    ClothMesh(const ObjMeshD& mesh,MESH_TYPE type);
    ClothMesh(const ClothMesh& other);
    ClothMesh& operator=(const ClothMesh& other);
	void reset(const ObjMeshD& mesh,MESH_TYPE type);
	bool write(std::ostream& os) const;
	bool read(std::istream& is);
	void writeVTKC(const std::string& str,char defColor=-1,const std::vector<char>* color=NULL,bool last=false) const;
    void writeVTKN(const std::string& str) const;
    void assignMass(MASS_MODE mode);
	void saveLast();
    void findBoundary(std::vector<std::vector<boost::shared_ptr<ClothTriangle> > >& boundary);
    void assembleIndex();
    void assembleA();
    void assembleN(Vec* N=NULL,Vec* NV=NULL,Vec* M=NULL,bool pos0=false) const;
    void assembleC(Vec* C=NULL,Vec* CV=NULL,Vec* M=NULL,bool pos0=false) const;
    void assignN(const Vec* N=NULL,const Vec* NV=NULL);
    void assignC(const Vec* C=NULL,const Vec* CV=NULL);
    void convertN2C();
    void convertC2N();
    void convertC2Obj(ObjMeshD& mesh);
    static Vec3d perp(const Vec3d& in);
    static void addI3x3(std::vector<Eigen::Triplet<scalarD,sizeType> >& trips,sizeType R,sizeType C,const scalarD& coef);
    template <typename T>
    static void add1x3(std::vector<Eigen::Triplet<scalarD,sizeType> >& trips,sizeType R,sizeType C,const T& v) {
        trips.push_back(Eigen::Triplet<scalarD,sizeType>(R,C+0,v[0]));
        trips.push_back(Eigen::Triplet<scalarD,sizeType>(R,C+1,v[1]));
        trips.push_back(Eigen::Triplet<scalarD,sizeType>(R,C+2,v[2]));
    }
    template <typename T>
    static void add3x3(std::vector<Eigen::Triplet<scalarD,sizeType> >& trips,sizeType R,sizeType C,const T& m) {
        for(sizeType r=0; r<3; r++)
            for(sizeType c=0; c<3; c++)
                trips.push_back(Eigen::Triplet<scalarD,sizeType>(R+r,C+c,m(r,c)));
    }
    void parityCheck();
    void parityCheckVss();
	//for rigid body
	Mat4d getTransform() const;
	void transform(const Mat3d& R,const Vec3d& T);
    //data
    std::vector<boost::shared_ptr<ClothVertex> > _vss;
    std::vector<boost::shared_ptr<ClothEdge> > _ess;
    std::vector<boost::shared_ptr<ClothTriangle> > _tss;
private:
    //fixed sparse matrix
    Eigen::SparseMatrix<scalarD,0,sizeType> _A,_AAT;
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<scalarD,0,sizeType> > _AATSol;
};

#endif
