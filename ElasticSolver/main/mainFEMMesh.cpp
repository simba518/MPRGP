#include "mainConfig.h"
#ifdef MAIN_FEMMESH

#include "FEMMesh.h"
#include "FEMGeom.h"
#include "FEMCollision.h"
#include "FEMCollider.h"
#include "ImplicitFuncInterface.h"
#include "CollisionDetection.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

class FuncSolid2D : public ImplicitFunc<scalar>
{
public:
    scalar operator()(const Vec3& pos) const {
        scalar minDist=numeric_limits<scalar>::max();
        for(scalar x=-9.0f; x<9.0f; x+=3.0f)
            for(scalar y=-9.0f; y<9.0f; y+=3.0f) {
                Vec3 posA=pos;
                posA.block<2,1>(0,0)+=Vec2(x,y);
                scalar dist=
                    std::max(std::max(0.1f-posA.x(),posA.x()-0.9f),
                             std::max(0.1f-posA.y(),posA.y()-0.9f));
                minDist=std::min(dist,minDist);
            }
        return minDist;
    }
};
class FuncSolid3D : public ImplicitFunc<scalar>
{
public:
    scalar operator()(const Vec3& pos) const {
        scalar minDist=numeric_limits<scalar>::max();
        for(scalar x=-9.0f; x<9.0f; x+=6.0f)
            for(scalar y=-9.0f; y<9.0f; y+=6.0f)
                for(scalar z=-9.0f; z<9.0f; z+=6.0f) {
                    Vec3 posA=pos+Vec3(x,y,z);
                    scalar dist=
                        std::max(std::max(std::max(0.1f-posA.x(),posA.x()-0.9f),
                                          std::max(0.1f-posA.y(),posA.y()-0.9f)),
                                 std::max(0.1f-posA.z(),posA.z()-0.9f));
                    minDist=std::min(dist,minDist);
                }
        return minDist;
    }
};
void mainFEMMesh1(std::string name="FEMMesh1")
{
    boost::filesystem::create_directory("./"+name);
    FEMMesh mesh(2,boost::shared_ptr<FEMCollision>(new FEMCollision));
    mesh.reset(BBox<scalar>(Vec3(-10.0f,-10.0f,0.0f),Vec3(10.0f,10.0f,0.0f)),FuncSolid2D(),0.06f,0.1f);
    mesh.writeVTK("./"+name+"/meshInit.vtk");
    mesh.writePSetVTK("./"+name+"/pSetInit.vtk");

    FEMMesh mesh2=mesh;
    Mat4 R=Mat4::Identity();
    R.block<3,1>(0,3)=Vec3(20.0f,20.0f,0.0f);
    mesh2.applyTrans(R,-1,true,true);
    mesh+=mesh2;
    mesh.writeVTK("./"+name+"/meshInit2.vtk");
    mesh.writePSetVTK("./"+name+"/pSetInit2.vtk");

    do {
        mesh-=rand()%mesh.nrB();
        ostringstream oss;
        oss << mesh.nrB();
        mesh.writeVTK("./"+name+"/mesh"+oss.str()+".vtk");
        mesh.writePSetVTK("./"+name+"/pSet"+oss.str()+".vtk");
    } while(mesh.nrB() > 1);
}
void mainFEMMesh2(std::string name="FEMMesh2")
{
    boost::filesystem::create_directory("./"+name);
    FEMMesh mesh(3,boost::shared_ptr<FEMCollision>(new FEMCollision));
    mesh.reset(BBox<scalar>(Vec3(-10.0f,-10.0f,-10.0f),Vec3(10.0f,10.0f,10.0f)),FuncSolid3D(),0.06f,0.1f);
    mesh.writeVTK("./"+name+"/meshInit.vtk");
    mesh.writePSetVTK("./"+name+"/pSetInit.vtk");

    FEMMesh mesh2=mesh;
    Mat4 R=Mat4::Identity();
    R.block<3,1>(0,3)=Vec3(20.0f,20.0f,20.0f);
    mesh2.applyTrans(R,-1,true,true);
    mesh+=mesh2;
    mesh.writeVTK("./"+name+"/meshInit2.vtk");
    mesh.writePSetVTK("./"+name+"/pSetInit2.vtk");

    do {
        mesh-=rand()%mesh.nrB();
        ostringstream oss;
        oss << mesh.nrB();
        mesh.writeVTK("./"+name+"/mesh"+oss.str()+".vtk");
        mesh.writePSetVTK("./"+name+"/pSet"+oss.str()+".vtk");
    } while(mesh.nrB() > 1);
}
void mainFEMMesh3(std::string name="FEMMesh3")
{
    boost::filesystem::create_directory("./"+name);
    FEMMesh mesh(3,boost::shared_ptr<FEMCollision>(new FEMCollision));
    mesh.reset("./mesh.ABQ",1.0f);
    mesh.writeVTK("./"+name+"/meshInit.vtk");
    mesh.writePSetVTK("./"+name+"/pSetInit.vtk");
}
void mainFEMMesh4(sizeType dim,boost::shared_ptr<FEMCollision> coll,std::string name="FEMMesh4")
{
    boost::filesystem::create_directory("./"+name);
    FEMMesh mesh(dim,coll);
    ParticleSetN pSet;
    scalar rad=0.1f;
    for(scalar k=(dim==2)?0.0f:-1.0f; k<1.0f; k+=(dim==2)?2.0f:rad)
        for(scalar i=-1.0f; i<1.0f; i+=rad)
            for(scalar j=-1.0f; j<1.0f; j+=rad) {
                ParticleN<scalar> p,pp;
                p._pos=Vec3(i,j,k);
                if(p._pos.norm() < 0.75f) {
                    for(scalar kk=dim==2.0f?0.0f:-1.0f; kk<=1.0f; kk+=2.0f)
                        for(scalar ii=-1.0f; ii<=1.0f; ii+=2.0f)
                            for(scalar jj=-1.0f; jj<=1.0f; jj+=2.0f) {
                                pp=p;
                                pp._pos[0]+=ii;
                                pp._pos[1]+=jj;
                                pp._pos[2]+=kk;
                                pSet.addParticle(pp);
                            }
                }
            }
    pSet.writeVTK("./"+name+"/pset.vtk");

    mesh.reset(0.2f,pSet);
    mesh.writeVTK("./"+name+"/mesh.vtk");

#define NR_TEST 100
    boost::filesystem::create_directory("./"+name+"/Raycast");
    for(sizeType i=0; i<NR_TEST; i++) {
        Vec3 r=Vec3::Zero();
        r.block(0,0,dim,1).setRandom();
#if DIM==2
        r[2]=0.0f;
#endif
        r=r.normalized()*4.0f;
        LineSeg l(r,Vec3::Zero());
        FEMInterp cd,cd2;
        scalar dist;

        {
            dist=numeric_limits<scalar>::max();
            mesh.getColl().rayCastMesh(l,cd,dist);
            if(dist < numeric_limits<scalar>::max())INFOV("%ld",i)

                vector<Vec3,Eigen::aligned_allocator<Vec3> > vss;
            vss.push_back(r);
            if(dist == numeric_limits<scalar>::max())
                vss.push_back(Vec3::Zero());
            else vss.push_back(cd._cell->getVert(cd));

            ostringstream oss;
            oss << "./"+name+"/Raycast/fmesh" << i << ".vtk";
            VTKWriter<scalar> os("Raycast",oss.str(),true);
            os.appendPoints(vss.begin(),vss.end());
            os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,2,0),
                           VTKWriter<scalar>::IteratorIndex<Vec3i>(1,2,0),
                           VTKWriter<scalar>::LINE);
        }

        {
            dist=numeric_limits<scalar>::max();
            mesh.getColl().rayCastPSet(l,rad,cd,dist);
            if(dist < numeric_limits<scalar>::max())INFOV("%ld",i)

                ostringstream oss;
            oss << "./"+name+"/Raycast/fpset" << i << ".vtk";
            VTKWriter<scalar> os("Raycast",oss.str(),true);
            if(dist < numeric_limits<scalar>::max()) {
                vector<Vec3,Eigen::aligned_allocator<Vec3> > vss;
                vss.push_back(cd._cell->getVert(cd));
                os.appendPoints(vss.begin(),vss.end());
                os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,0,0),
                               VTKWriter<scalar>::IteratorIndex<Vec3i>(1,0,0),
                               VTKWriter<scalar>::POINT);
            }
        }
    }
}
void mainFEMMesh5(boost::shared_ptr<FEMCollision> coll,std::string name="FEMMesh5")
{
    boost::filesystem::create_directory("./"+name);
    FEMGeom geom(3);
    Mat4 R=Mat4::Identity();
    for(sizeType i=0; i<5; i++) {
        geom.addGeomMesh(R,"./bunny.obj",0.1f);
        R.block<3,1>(0,3).array()+=3.0f;
    }
    geom.assemble();
    geom.writeVTK("./"+name+"/CGeometry3Dgeom.vtk");

#define SZ 0.1f
#define RAD 0.5f
    ParticleSetN pSet;
    for(scalar k=-RAD; k<RAD; k+=SZ)
        for(scalar i=-RAD; i<RAD; i+=SZ)
            for(scalar j=-RAD; j<RAD; j+=SZ) {
                ParticleN<scalar> p;
                p._pos=Vec3(i,j,k);
                if(p._pos.norm() < RAD) {
                    p._pos.array()-=RAD*5.0f;
                    pSet.addParticle(p);
                }
            }
    pSet.writeVTK("./"+name+"/CGeometry3Dpset.vtk");
    FEMMesh mesh(3,coll);
    mesh.reset(0.2f,pSet);
    mesh.writeVTK("./"+name+"/CGeometry3Dmesh.vtk");

    boost::filesystem::create_directory("./"+name+"/CollGeomCD");
    for(scalar i=0,j=0; i<20.0f; i+=0.25f) {
        Vec3 dv=Vec3::Zero();
        dv.setConstant(i);
        for(sizeType b=0; b<mesh.nrB(); b++)
            for(sizeType c=0; c<mesh.getB(b).nrV(); c++) {
                FEMVertex& v=mesh.getB(b).getV(c);
                v._pos=v._pos0+dv;
            }
        mesh.updateMesh();

        ostringstream oss;
        oss << "./"+name+"/CollGeomCD/mesh" << j << ".vtk";
        mesh.writeVTK(oss.str());

        ostringstream oss2;
        oss2 << "./"+name+"/CollGeomCD/pSet" << j << ".vtk";
        mesh.writePSetVTK(oss2.str());

        ostringstream oss3;
        oss3 << "./"+name+"/CollGeomCD/coll" << j << ".vtk";
        DebugFEMCollider coll(oss3.str(),3);
        mesh.getColl().collideGeom(geom,coll,true);
        j++;
    }
}
void mainFEMMesh6(sizeType dim,boost::shared_ptr<FEMCollision> coll,std::string name="FEMMesh6")
{
    boost::filesystem::create_directory("./"+name);
    FEMMesh mesh(dim,coll);
#define rad 0.05f
    ParticleSetN pSet;
    for(scalar k=(dim==2)?0.0f:-1.0f; k<1.0f; k+=(dim==2)?2.0f:rad)
        for(scalar i=-1.0f; i<1.0f; i+=rad)
            for(scalar j=-1.0f; j<1.0f; j+=rad) {
                Vec3 dv=Vec3::Zero();
                dv.block(0,0,dim,1).setConstant(3.0f);

                ParticleN<scalar> p;
                p._pos=Vec3(i,j,k);
                pSet.addParticle(p);
                p._pos+=dv;
                pSet.addParticle(p);
            }
    pSet.writeVTK("./"+name+"/pset.vtk");

    srand(1000);
    mesh.reset(0.2f,pSet);
    Mat4 T=Mat4::Identity();
    T.block<3,3>(0,0)=makeRotation<scalar>(Vec3(0.0f,0.0f,(scalar)rand()));
    mesh.applyTrans(T,0,true,true);
    mesh.updateMesh();
    mesh.writeVTK("./"+name+"/mesh.vtk");

    boost::filesystem::create_directory("./"+name+"/CollSelfCD/");
    sizeType j=0;
    for(scalar i=0.0f; i<5.0f; i+=0.02f) {
        Vec3 dv=Vec3::Zero();
        dv.block(0,0,dim,1).setConstant(i);
        for(sizeType v=0; v<(sizeType)mesh.getB(0).nrV(); v++)
            mesh.getB(0).getV(v)._pos=mesh.getB(0).getV(v)._pos0+dv;
        mesh.updateMesh();

        ostringstream oss;
        oss << "./"+name+"/CollSelfCD/mesh" << j << ".vtk";
        mesh.writeVTK(oss.str());

        ostringstream oss2;
        oss2 << "./"+name+"/CollSelfCD/pSet" << j << ".vtk";
        mesh.writePSetVTK(oss2.str());

        ostringstream oss3;
        oss3 << "./"+name+"/CollSelfCD/coll" << j << ".vtk";
        DebugFEMCollider coll(oss3.str(),dim);
        mesh.getColl().collideMesh(coll,false);

        ostringstream oss4;
        oss4 << "./"+name+"/CollSelfCD/collS" << j << ".vtk";
        DebugFEMCollider collS(oss4.str(),dim);
        mesh.getColl().collideMesh(collS,true);
        j++;
    }
}
void main()
{
    //mainFEMMesh1();
    //mainFEMMesh2();
    //mainFEMMesh3();
    //mainFEMMesh4(2,boost::shared_ptr<FEMCollision>(new FEMCollision),"FEMRayCast2D");
    //mainFEMMesh4(3,boost::shared_ptr<FEMCollision>(new FEMCollision),"FEMRayCast3D");
    //mainFEMMesh4(2,boost::shared_ptr<FEMCollision>(new BVHFEMCollision),"FEMRayCast2DBVH");
    //mainFEMMesh4(3,boost::shared_ptr<FEMCollision>(new BVHFEMCollision),"FEMRayCast3DBVH");
    //mainFEMMesh4(2,boost::shared_ptr<FEMCollision>(new SBVHFEMCollision),"FEMRayCast2DSBVH");
    //mainFEMMesh4(3,boost::shared_ptr<FEMCollision>(new SBVHFEMCollision),"FEMRayCast3DSBVH");
    //mainFEMMesh5(boost::shared_ptr<FEMCollision>(new FEMCollision),"FEMCGeom3D");
    //mainFEMMesh5(boost::shared_ptr<FEMCollision>(new BVHFEMCollision),"FEMCGeom3DBVH");
    //mainFEMMesh5(boost::shared_ptr<FEMCollision>(new SBVHFEMCollision),"FEMCGeom3DSBVH");
    //mainFEMMesh6(2,boost::shared_ptr<FEMCollision>(new FEMCollision),"FEMCollideMesh2D");
    //mainFEMMesh6(2,boost::shared_ptr<FEMCollision>(new BVHFEMCollision),"FEMCollideMesh2DBVH");
    //mainFEMMesh6(2,boost::shared_ptr<FEMCollision>(new SBVHFEMCollision),"FEMCollideMesh2DSBVH");
    //mainFEMMesh6(3,boost::shared_ptr<FEMCollision>(new FEMCollision),"FEMCollideMesh3D");
    //mainFEMMesh6(3,boost::shared_ptr<FEMCollision>(new BVHFEMCollision),"FEMCollideMesh3DBVH");
    //mainFEMMesh6(3,boost::shared_ptr<FEMCollision>(new SBVHFEMCollision),"FEMCollideMesh3DSBVH");
}

#endif
