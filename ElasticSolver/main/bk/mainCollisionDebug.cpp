#include "FEMSolver.h"
#include "FEMMesh.h"
#include "FEMGeom.h"
#include "FEMCollision.h"
#include "FEMCollider.h"
#include "FEMProfiler.h"
#include "FEMLocalBasis.h"
#include "MakeMesh.h"
#include "CollisionDetection.h"
#include "ImplicitFuncInterface.h"
#include <boost/shared_ptr.hpp>
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE
	
#define DIM 3
#define SZ 0.1f
#define BODYID 0
void mainRaycast(boost::shared_ptr<FEMCollision> coll)
{
//#define SINGLE_BALL
	ParticleSetN pSet;
#if DIM==2
	scalar k=0.0f;
#else
	for(scalar k=-1;k<1;k+=SZ)
#endif
	for(scalar i=-1;i<1;i+=SZ)
	for(scalar j=-1;j<1;j+=SZ)
	{
		ParticleN<scalar> p,pp;
		p._pos=Vec3(i,j,k);
#ifdef SINGLE_BALL
		if(p._pos.norm() < 1)
			pSet.addParticle(p);
#else	
		if(p._pos.norm() < 0.75)
		{
#if DIM==2
			scalar kk=0.0f;
#else
			for(scalar kk=-1;kk<=1;kk+=2)
#endif
			for(scalar ii=-1;ii<=1;ii+=2)
			for(scalar jj=-1;jj<=1;jj+=2)
			{
				pp=p;
				pp._pos[0]+=ii;
				pp._pos[1]+=jj;
				pp._pos[2]+=kk;
				pSet.addParticle(pp);
			}
		}
#endif
	}
	pSet.writeVTK("./pset.vtk");

	FEMMesh mesh(DIM,coll);
	mesh.reset(DIM,0.2f,pSet);
	mesh.writeVTK("./mesh.vtk");

#define NR_TEST 100
	boost::filesystem::create_directory("./Raycast");
	for(sizeType i=0;i<NR_TEST;i++)
	{
		Vec3 r=Vec3::Random();
#if DIM==2
		r[2]=0.0f;
#endif
		r=r.normalized()*3;
		LineSeg l(r,Vec3::Zero());
		FEMInterp cd,cd2;
		scalar dist;

		{
			dist=numeric_limits<scalar>::max();
			mesh.getColl().rayCastMesh(l,cd,dist,BODYID);
			if(dist < numeric_limits<scalar>::max())INFOV("%ld",i)

			vector<Vec3,Eigen::aligned_allocator<Vec3> > vss;
			vss.push_back(r);
			if(dist == numeric_limits<scalar>::max())
				vss.push_back(Vec3::Zero());
			else vss.push_back(mesh.getC(cd._id).getVert(cd));

			ostringstream oss;oss << "./Raycast/fmesh" << i << ".vtk";
			VTKWriter<scalar> os("Raycast",oss.str(),true);
			os.appendPoints(vss.begin(),vss.end());
			os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,2,0),
						   VTKWriter<scalar>::IteratorIndex<Vec3i>(1,2,0),
						   VTKWriter<scalar>::LINE);
		}

		{
			dist=numeric_limits<scalar>::max();
			mesh.getColl().rayCastPSet(l,SZ,cd,dist,BODYID);
			if(dist < numeric_limits<scalar>::max())INFOV("%ld",i)

			ostringstream oss;oss << "./Raycast/fpset" << i << ".vtk";
			VTKWriter<scalar> os("Raycast",oss.str(),true);
			if(dist < numeric_limits<scalar>::max())
			{
				vector<Vec3,Eigen::aligned_allocator<Vec3> > vss;
				vss.push_back(mesh.getPSet()[cd._id]._pos);
				os.appendPoints(vss.begin(),vss.end());
				os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,0,0),
							   VTKWriter<scalar>::IteratorIndex<Vec3i>(1,0,0),
							   VTKWriter<scalar>::POINT);
			}
		}
	}
}
void mainCGeometry3D(boost::shared_ptr<FEMCollision> coll)
{
	FEMGeom geom(3);
	geom.addGeomMesh(Mat4::Identity(),"./bunny.obj",0.1f);
	geom.assemble();
	geom.writeVTK("./CGeometry3Dgeom.vtk");

#define SZ 0.1f
#define RAD 0.5f
	ParticleSetN pSet;
	for(scalar k=-RAD;k<RAD;k+=SZ)
	for(scalar i=-RAD;i<RAD;i+=SZ)
	for(scalar j=-RAD;j<RAD;j+=SZ)
	{
		ParticleN<scalar> p;
		p._pos=Vec3(i,j,k);
		if(p._pos.norm() < RAD)
		{
			p._pos.array()-=RAD*5.0f;
			pSet.addParticle(p);
		}
	}
	pSet.writeVTK("./CGeometry3Dpset.vtk");
	FEMMesh mesh(3,coll);
	mesh.reset(3,0.2f,pSet);
	mesh.writeVTK("./CGeometry3Dmesh.vtk");

	boost::filesystem::create_directory("./CollGeomCD");
	for(scalar i=0,j=0;i<RAD*5.0f;i+=0.01f)
	{
		Vec3 dv=Vec3::Zero();
		dv.setConstant(i);
		for(sizeType c=0;c<(sizeType)mesh.nrV();c++){
			FEMVertex& v=mesh.getV(c);
			v._pos=v._pos0+dv;
		}
		mesh.updateMesh();
	
		ostringstream oss;oss << "./CollGeomCD/mesh" << j << ".vtk";
		mesh.writeVTK(oss.str());

		ostringstream oss2;oss2 << "./CollGeomCD/pSet" << j << ".vtk";
		mesh.getPSet().writeVTK(oss2.str());
		
		ostringstream oss3;oss3 << "./CollGeomCD/coll" << j << ".vtk";
		DebugFEMCollider coll(oss3.str(),3);
		mesh.getColl().collideGeom(geom,coll,true);
		j++;
	}
}
void mainGeometry2D()
{
	FEMGeom geom(2);
	//box
	scalar theta=rand()*M_PI/(scalar)RAND_MAX;
	Mat2 r;r << cos(theta),sin(theta),-sin(theta),cos(theta);
	geom.addGeomBox(OBB2D(r,Vec2::Random(),Vec2::Random()+Vec2::Ones()));
	//plane
	Vec3 x0=Vec3::Random();x0[2]=0.0f;
	Vec3 n=Vec3::Random();n[2]=0.0f;
	geom.addGeomPlane(Mat4::Identity(),Plane(x0,n));
	//sphere
	for(scalar i=10;i<30;i+=5)
	for(scalar j=10;j<30;j+=5)
		geom.addGeomSphere(Vec3(i,j,0.0f),1.5f);
	//mesh
	ObjMesh mesh;
	MakeMesh::makeCapsule2D(mesh,1.0f,5.0f,16);
	geom.addGeomMesh(Mat4::Identity(),mesh,0.5f);
	geom.assemble();
	//write
	geom.write(boost::filesystem::ofstream("./Geometry2Dgeom.dat",ios::binary));
	
	FEMGeom geom2(2);
    geom2.read(boost::filesystem::ifstream("./Geometry2Dgeom.dat",ios::binary));
	geom2.writeVTK("./Geometry2Dgeom.vtk");
	geom2.writeBVH();
}
void mainGeometry3D()
{
	FEMGeom geom(3);
	//box
	geom.addGeomBox(OBB3D(makeRotation<scalar>(Vec3::Random()),
						  Vec3::Random(),
						  BBox<scalar>(Vec3::Random()-Vec3::Ones(),
									   Vec3::Random()+Vec3::Ones())));
	//plane
	geom.addGeomPlane(Mat4::Identity(),Vec4(0,1,1,0.0f));
	geom.addGeomPlane(Mat4::Identity(),Plane(Vec3(0,50,50),Vec3(0,-1,-1)));
	//sphere
	for(scalar i=10;i<30;i+=5)
	for(scalar j=10;j<30;j+=5)
	for(scalar k=10;k<30;k+=5)
		geom.addGeomSphere(Vec3(i,j,k),1.5f);
	geom.addGeomMesh(Mat4::Identity(),"./bunny.obj",0.1f);
	//mesh
	ObjMesh mesh;
	MakeMesh::makeCapsule3D(mesh,1.0f,5.0f,16);
	geom.addGeomMesh(Mat4::Identity(),mesh,0.5f);
	geom.assemble();
	//write
	geom.write(boost::filesystem::ofstream("./Geometry3Dgeom.dat",ios::binary));
	
	FEMGeom geom2(3);
	geom2.read(boost::filesystem::ifstream("./Geometry3Dgeom.dat",ios::binary));
	geom2.writeVTK("./Geometry3Dgeom.vtk");
	geom2.writeBVH();
}
void mainCollideMesh(boost::shared_ptr<FEMCollision> coll)
{
#define RAD 0.5f
	ParticleSetN pSet;
#if DIM==2
	scalar k=0.0f;
#else
	for(scalar k=-RAD;k<RAD;k+=SZ)
#endif
	for(scalar i=-RAD;i<RAD;i+=SZ)
	for(scalar j=-RAD;j<RAD;j+=SZ)
	{
		Vec3 dv=Vec3::Zero();
		dv.block(0,0,DIM,1).setConstant(RAD);

		ParticleN<scalar> p;
		p._pos=Vec3(i,j,k);
		if(p._pos.norm() < RAD)
		{
			pSet.addParticle(p);
			p._pos+=dv*2.5f;
			pSet.addParticle(p);
		}
	}
	pSet.writeVTK("./pset.vtk");

	FEMMesh mesh(DIM,coll);
	mesh.reset(DIM,0.2f,pSet);
	Mat4 T=Mat4::Identity();
	T.block<3,3>(0,0)=makeRotation<scalar>(Vec3(0.0f,0.0f,(scalar)rand()));
	mesh.applyTrans(T,0,true,true);
	mesh.updateMesh();
	mesh.writeVTK("./mesh.vtk");

	boost::filesystem::create_directory("./CollSelfCD/");
	sizeType j=0;
	for(scalar i=0;i<RAD*5.0f;i+=0.02f)
	{
		Vec3 dv=Vec3::Zero();
		dv.block(0,0,DIM,1).setConstant(i);
		for(sizeType c=0;c<(sizeType)mesh.nrC();c++)
			if(mesh.getC(c)._bodyId==BODYID)
			for(sizeType vid=0;vid<DIM+1;vid++)
			{
				boost::shared_ptr<FEMVertex> v=mesh.getC(c)._v[vid];
				v->_pos=v->_pos0+dv;
			}
		mesh.updateMesh();
		
		ostringstream oss;oss << "./CollSelfCD/mesh" << j << ".vtk";
		mesh.writeVTK(oss.str());

		ostringstream oss2;oss2 << "./CollSelfCD/pSet" << j << ".vtk";
		mesh.getPSet().writeVTK(oss2.str());
		
		ostringstream oss3;oss3 << "./CollSelfCD/coll" << j << ".vtk";
		DebugFEMCollider coll(oss3.str(),DIM);
		mesh.getColl().collideMesh(coll,false);

		ostringstream oss4;oss4 << "./CollSelfCD/collS" << j << ".vtk";
		DebugFEMCollider collS(oss4.str(),DIM);
		mesh.getColl().collideMesh(collS,true);
		j++;
	}
}
void mainBodyOp(boost::shared_ptr<FEMCollision> coll)
{
#define RAD 0.5f
	ParticleSetN pSetBall;
#if DIM==2
	scalar k=0.0f;
#else
	for(scalar k=-RAD;k<RAD;k+=RAD/2)
#endif
	for(scalar i=-RAD;i<RAD;i+=RAD/2)
	for(scalar j=-RAD;j<RAD;j+=RAD/2)
	{
		ParticleN<scalar> p;
		p._pos=Vec3(i,j,k);
		if(p._pos.norm() < RAD)
			pSetBall.addParticle(p);
	}
#define SZS 2.0f
#define SP 5.0f
	ParticleSetN pSet;
#if DIM==2
#else
	for(scalar k=0;k<SP;k+=SZS)
#endif
	for(scalar i=0;i<SP;i+=SZS)
	for(scalar j=0;j<SP;j+=SZS)
	{
		Vec3 rand=Vec3::Zero();
		rand.block<2,1>(0,0).setRandom();
		rand*=0.25f;
		for(sizeType p=0;p<pSetBall.size();p++)
		{
			ParticleN<scalar> pp=pSetBall[p];
			pp._pos+=Vec3(i,j,k)+rand;
			pSet.addParticle(pp);
		}
	}
	pSet.writeVTK("./psetBodyOp.vtk");

	{
		FEMMesh mesh(DIM,coll);
		mesh.reset(DIM,0.4f,pSet);
		mesh.writeVTK("./meshBodyOp.vtk");

		mesh-=rand()%mesh.nrB();
		mesh.writeVTK("./meshBodyOp_1.vtk");
		mesh-=rand()%mesh.nrB();
		mesh.writeVTK("./meshBodyOp_2.vtk");

		FEMMesh meshTmp=mesh;
		Mat4 R=Mat4::Identity();
		R.block<3,1>(0,3)=Vec3::Constant(100.0f);
		meshTmp.applyTrans(R,-1);
		mesh+=meshTmp;
		mesh.write(boost::filesystem::ofstream("./tmp.dat",ios::binary));
	}

	FEMMesh mesh2(DIM,coll->copy());
	mesh2.read(boost::filesystem::ifstream("./tmp.dat",ios::binary));
	mesh2.writeVTK("./meshBodyOp_3.vtk");

	boost::filesystem::create_directory("./Raycast2");
	for(sizeType i=0;i<NR_TEST;i++)
	{
		Vec3 r=Vec3::Random();
#if DIM==2
		r[2]=0.0f;
#endif
		r=r.normalized()*5;
		Vec3 r0=Vec3::Constant(2.5f);
		LineSeg l(r+r0,r0);
		FEMInterp cd,cd2;
		scalar dist;

		{
			dist=numeric_limits<scalar>::max();
			mesh2.getColl().rayCastMesh(l,cd,dist,-1);
			if(dist < numeric_limits<scalar>::max())INFOV("%ld",i)

			vector<Vec3,Eigen::aligned_allocator<Vec3> > vss;
			vss.push_back(r+r0);
			if(dist == numeric_limits<scalar>::max())
				vss.push_back(r0);
			else vss.push_back(mesh2.getC(cd._id).getVert(cd));

			ostringstream oss;oss << "./Raycast2/fmesh" << i << ".vtk";
			VTKWriter<scalar> os("Raycast2",oss.str(),true);
			os.appendPoints(vss.begin(),vss.end());
			os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,2,0),
						   VTKWriter<scalar>::IteratorIndex<Vec3i>(1,2,0),
						   VTKWriter<scalar>::LINE);
		}

		{
			dist=numeric_limits<scalar>::max();
			mesh2.getColl().rayCastPSet(l,SZ,cd,dist,-1);
			if(dist < numeric_limits<scalar>::max())INFOV("%ld",i)

			ostringstream oss;oss << "./Raycast2/fpset" << i << ".vtk";
			VTKWriter<scalar> os("Raycast2",oss.str(),true);
			if(dist < numeric_limits<scalar>::max())
			{
				vector<Vec3,Eigen::aligned_allocator<Vec3> > vss;
				vss.push_back(mesh2.getPSet()[cd._id]._pos);
				os.appendPoints(vss.begin(),vss.end());
				os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,0,0),
							   VTKWriter<scalar>::IteratorIndex<Vec3i>(1,0,0),
							   VTKWriter<scalar>::POINT);
			}
		}
	}
}

class FuncSolid : public ImplicitFunc<scalar>
{
public:
	FuncSolid(sizeType dim):_dim(dim){}
	scalar operator()(const Vec3& pos) const
	{
		return std::max(std::max(0.1f-pos.x(),pos.x()-9.3f),
						std::max(0.1f-pos.y(),pos.y()-9.3f));
	}
	sizeType _dim;
};
int main()
{
	boost::shared_ptr<FEMCollision> coll(new FEMCollision);

	//FEMMesh mesh(2,coll);
	//FuncSolid f(2);
	//mesh.reset(BBox<scalar>(Vec3::Zero(),Vec3(10.0f,10.0f,0.0f)),f,0.1f,0.2f);

	FEMMesh mesh(3,coll);
//#define READ
#ifndef READ
	mesh.reset(boost::filesystem::ifstream("F:/data/cheb/mesh.ABQ"),0.0f);
	mesh.write(boost::filesystem::ofstream("F:/data/cheb/mesh.dat",ios::binary));
#else
	mesh.read(boost::filesystem::ifstream("F:/data/cheb/mesh.dat",ios::binary));
#endif

	mesh.getPSet().writeVTK("./pset.vtk");
	mesh.writeVTK("./mesh.vtk");
	FEMMesh mesh2(3,coll->copy());
	mesh2=mesh;
	Mat4 T=Mat4::Identity();
	T.block<2,1>(0,3).setConstant(100.0f);
	mesh2.applyTrans(T,-1);
	mesh+=mesh2;

	FEMLocalBasis basis;
	basis.updateMesh(mesh.getB(1));
	//basis.writeVTK(mesh.getB(1),"./tyran/");
	basis.debugPatchVTK(mesh.getB(1),false,10);
}
/*int main()
{
	boost::shared_ptr<FEMCollision> coll(new SBVHFEMCollision);
	//mainGeometry2D();
	//mainGeometry3D();
	//mainCGeometry3D(coll);
	//mainRaycast(coll);
	mainCollideMesh(coll);
	//mainBodyOp(coll);
}*/