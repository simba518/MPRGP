#ifndef CLOTH_COLLISION_H
#define CLOTH_COLLISION_H

#include <geom/BVHBuilder.h>
#include <CollisionDetection.h>
#include "ClothMesh.h"

USE_PRJ_NAMESPACE

class ClothCollision
{
public:
	//typedef BBox<scalar> BBOX;
	typedef KDOP18<scalar> BBOX;
	struct CollisionHandler 
	{
		CollisionHandler();
		virtual ~CollisionHandler(){}
		virtual void handle(boost::shared_ptr<ClothMesh::ClothVertex> V1,boost::shared_ptr<ClothMesh::ClothTriangle> T2,const Vec3d n,const Vec4d& omg,scalarD t);
        virtual void handle(boost::shared_ptr<ClothMesh::ClothEdge> E1,boost::shared_ptr<ClothMesh::ClothEdge> E2,const Vec3d n,const Vec4d& omg,scalarD t);
	private:
		sizeType _id;
	};
	struct NarrowNode : public Serializable
	{
		NarrowNode();
		void buildBVH();
		BBOX refit(bool useLastPos=true);
		bool read(std::istream& is,IOData* dat);
		bool write(std::ostream& os,IOData* dat) const;
		bool operator<(const NarrowNode& other) const;
		vector<BBOX> _vbb,_ebb;
		boost::shared_ptr<ClothMesh> _mesh;
	    vector<GEOM::Node<sizeType,BBOX> > _bvh;
		vector<bool> _active;
	};
	template <typename TA,typename TB>
	struct Cache
	{
		Cache(){}
		Cache(boost::shared_ptr<TA> A,boost::shared_ptr<TB> B):_A(A),_B(B){}
		bool operator==(const Cache& other) const{
			return _A==other._A && _B==other._B;
		}
		bool operator<(const Cache& other) const
		{
			if(_A<other._A)return true;
			else if(other._A<_A)return false;
			if(_B<other._B)return true;
			else if(other._B<_B)return false;
			return false;
		}
		boost::shared_ptr<TA> _A;
		boost::shared_ptr<TB> _B;
	};
	template <typename T>
	struct Cache<T,T>
	{
		Cache(){}
		Cache(boost::shared_ptr<T> A,boost::shared_ptr<T> B):_A(A),_B(B){if(_B<_A)std::swap(_A,_B);}
		bool operator==(const Cache& other) const{
			return _A==other._A && _B==other._B;
		}
		bool operator<(const Cache& other) const
		{
			if(_A<other._A)return true;
			else if(other._A<_A)return false;
			if(_B<other._B)return true;
			else if(other._B<_B)return false;
			return false;
		}
		boost::shared_ptr<T> _A;
		boost::shared_ptr<T> _B;
	};
public:
	sizeType nrMesh() const;
	const ClothMesh& getMesh(sizeType i) const;
	void updateMesh(boost::shared_ptr<ClothMesh> mesh);
	void addMesh(boost::shared_ptr<ClothMesh> mesh);
	void delMesh(boost::shared_ptr<ClothMesh> mesh);
	void collide(CollisionHandler& handler,bool useActive=false);
	void restartActive();
	void activate(const ClothMesh::ClothVertex* t);
	void onCell(const GEOM::Node<boost::shared_ptr<NarrowNode>,BBOX>& A,const GEOM::Node<boost::shared_ptr<NarrowNode>,BBOX>& B);
	void onCell(const GEOM::Node<sizeType,BBOX>& nA,const GEOM::Node<sizeType,BBOX>& nB);
private:
	//data
	vector<GEOM::Node<boost::shared_ptr<NarrowNode>,BBOX> > _bvh;
	vector<Cache<ClothMesh::ClothVertex,ClothMesh::ClothTriangle> > _cacheVT;
	vector<Cache<ClothMesh::ClothEdge,ClothMesh::ClothEdge> > _cacheEE;
	vector<Cache<NarrowNode,NarrowNode> > _cache;
	Cache<NarrowNode,NarrowNode> _activeCache;
public:
	//param
	static scalarD _thickness;
	static scalarD _rounding;
	static scalarD _timeRes;
};

#endif
