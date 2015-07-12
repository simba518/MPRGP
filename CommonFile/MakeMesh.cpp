#include "MakeMesh.h"

USE_PRJ_NAMESPACE

void MakeMesh::makeTet3D(ObjMesh& m){
    m.getV().clear();
    
    m.getV().push_back(Vec3(0.0f,0.0f,0.0f));
    m.getV().push_back(Vec3(1.0f,0.0f,0.0f));
    m.getV().push_back(Vec3(0.0f,1.0f,0.0f));
    m.getV().push_back(Vec3(0.0f,0.0f,1.0f));

    m.getI().push_back(Vec3i(0,2,1));
    m.getI().push_back(Vec3i(0,1,3));
    m.getI().push_back(Vec3i(0,3,2));
    m.getI().push_back(Vec3i(1,2,3));
}
void MakeMesh::makeBox3D(ObjMesh& m,const Vec3& ext){
    m.getV().clear();
#define ADDV(S1,S2,S3) m.getV().push_back(Vec3(S1*ext[0],S2*ext[1],S3*ext[2]));
    ADDV(-1.0f,-1.0f,-1.0f)
    ADDV(+1.0f,-1.0f,-1.0f)
    ADDV(+1.0f,+1.0f,-1.0f)
    ADDV(-1.0f,+1.0f,-1.0f)
    ADDV(-1.0f,-1.0f,+1.0f)
    ADDV(+1.0f,-1.0f,+1.0f)
    ADDV(+1.0f,+1.0f,+1.0f)
    ADDV(-1.0f,+1.0f,+1.0f)
    m.getI().clear();
#define ADDI(A,B,C) m.getI().push_back(Vec3i(A,B,C));
    //-Z
    ADDI(0,2,1)
    ADDI(0,3,2)
    //+Z
    ADDI(4,5,6)
    ADDI(4,6,7)
    //-X
    ADDI(0,4,7)
    ADDI(0,7,3)
    //-X
    ADDI(1,2,6)
    ADDI(1,6,5)
    //-Y
    ADDI(0,1,5)
    ADDI(0,5,4)
    //+Y
    ADDI(2,3,7)
    ADDI(2,7,6)
    m.applyTrans();
}
void MakeMesh::makeDiscreteBox3D(ObjMesh& m,const Vec3& ext){
    m.getV().clear();
	m.getI().clear();
	sizeType id;
#define ADDV(S1,S2,S3) m.getV().push_back(Vec3(S1*ext[0],S2*ext[1],S3*ext[2]));
	//top
	id=(sizeType)m.getV().size();
	ADDV(-1.0f,-1.0f,+1.0f)
    ADDV(+1.0f,-1.0f,+1.0f)
    ADDV(+1.0f,+1.0f,+1.0f)
    ADDV(-1.0f,+1.0f,+1.0f)
	ADDI(id+0,id+1,id+2)
	ADDI(id+0,id+2,id+3)
	//bottom
	id=(sizeType)m.getV().size();
	ADDV(-1.0f,-1.0f,-1.0f)
    ADDV(-1.0f,+1.0f,-1.0f)
    ADDV(+1.0f,+1.0f,-1.0f)
    ADDV(+1.0f,-1.0f,-1.0f)
	ADDI(id+0,id+1,id+2)
	ADDI(id+0,id+2,id+3)
	//right
	id=(sizeType)m.getV().size();
	ADDV(+1.0f,-1.0f,-1.0f)
    ADDV(+1.0f,+1.0f,-1.0f)
    ADDV(+1.0f,+1.0f,+1.0f)
    ADDV(+1.0f,-1.0f,+1.0f)
	ADDI(id+0,id+1,id+2)
	ADDI(id+0,id+2,id+3)
	//left
	id=(sizeType)m.getV().size();
	ADDV(-1.0f,-1.0f,-1.0f)
    ADDV(-1.0f,-1.0f,+1.0f)
    ADDV(-1.0f,+1.0f,+1.0f)
    ADDV(-1.0f,+1.0f,-1.0f)
	ADDI(id+0,id+1,id+2)
	ADDI(id+0,id+2,id+3)
	//front
	id=(sizeType)m.getV().size();
	ADDV(-1.0f,+1.0f,-1.0f)
    ADDV(-1.0f,+1.0f,+1.0f)
    ADDV(+1.0f,+1.0f,+1.0f)
    ADDV(+1.0f,+1.0f,-1.0f)
	ADDI(id+0,id+1,id+2)
	ADDI(id+0,id+2,id+3)
	//back
	id=(sizeType)m.getV().size();
	ADDV(-1.0f,-1.0f,-1.0f)
    ADDV(+1.0f,-1.0f,-1.0f)
    ADDV(+1.0f,-1.0f,+1.0f)
    ADDV(-1.0f,-1.0f,+1.0f)
	ADDI(id+0,id+1,id+2)
	ADDI(id+0,id+2,id+3)
	m.applyTrans();
}
void MakeMesh::makeBox2D(ObjMesh& m,const Vec3& ext){
    m.getV().clear();
    ADDV(-1.0f,-1.0f,0.0f)
    ADDV(+1.0f,-1.0f,0.0f)
    ADDV(+1.0f,+1.0f,0.0f)
    ADDV(-1.0f,+1.0f,0.0f)
    m.getI().clear();
    ADDI(0,1,0)
    ADDI(1,2,0)
    ADDI(2,3,0)
    ADDI(3,0,0)
    m.getDim()=2;
}
void MakeMesh::makeSphere3D(ObjMesh& m,const scalar& rad,const sizeType& slice){
    makeBox3D(m,Vec3::Ones());
    sizeType nr=(sizeType)std::ceil(log((scalar)slice));
    m.subdivide((int)nr);
    for(sizeType i=0;i<(sizeType)m.getV().size();i++)
        m.getV()[i]=m.getV()[i].normalized()*rad;
    m.applyTrans();
}
void MakeMesh::makeSphere2D(ObjMesh& m,const scalar& rad,const sizeType& slice){
    m.getV().clear();
    for(sizeType i=0;i<slice;i++){
        scalar ang=2.0f*M_PI*(i/(scalar)slice);
        m.getV().push_back(Vec3(cos(ang)*rad,sin(ang)*rad,0.0f));
    }
    m.getI().clear();
    for(sizeType i=0;i<slice;i++)
        m.getI().push_back(Vec3i(i,(i+1)%slice,0));
    m.getDim()=2;
}
void MakeMesh::makeCapsule3D(ObjMesh& m,const scalar& rad,const scalar& y,const sizeType& slice){
    makeSphere3D(m,rad,slice);
    for(sizeType i=0;i<(sizeType)m.getV().size();i++){
        Vec3& v=m.getV()[i];
        if(v.y() < 0.0f)
            v+=Vec3(0,-y,0);
        else v+=Vec3(0,y,0);
    }
    m.applyTrans();
}
void MakeMesh::makeCapsule2D(ObjMesh& m,const scalar& rad,const scalar& y,const sizeType& slice){
    makeSphere2D(m,rad,slice);
    for(sizeType i=0;i<(sizeType)m.getV().size();i++){
        Vec3& v=m.getV()[i];
        if(v.y() < 0.0f)
            v+=Vec3(0,-y,0);
        else v+=Vec3(0,y,0);
    }
    m.getDim()=2;
}
void MakeMesh::makeCylinder3D(ObjMesh& m,const scalar& rad,const scalar& y,const sizeType& slice,const sizeType& sliceY,bool cap){
    m.getV().clear();
	m.getI().clear();
	const scalar deltaY=y*2.0f/(scalar)sliceY;
	scalar Y;
	for(Y=-y;Y<=y+EPS;Y+=deltaY){
		sizeType IBeg=(sizeType)m.getV().size();
		for(sizeType i=0;i<slice;i++){
			if(IBeg > 0){
				sizeType I=IBeg+i;
				sizeType NI=IBeg+(i+1)%slice;
				m.getI().push_back(Vec3i(I,NI-slice,I-slice));
				m.getI().push_back(Vec3i(I,NI,NI-slice));
			}
		    scalar ang=2.0f*M_PI*((scalar)i/(scalar)slice);
		    m.getV().push_back(Vec3((scalar)cos(ang)*rad,Y,(scalar)sin(ang)*rad));
		}
	}
	if(cap){
		sizeType base=m.getV().size()-slice;
		m.getV().push_back(Vec3(0.0f,-y,0.0f));
		m.getV().push_back(Vec3(0.0f,Y-deltaY,0.0f));
		for(sizeType i=0;i<slice;i++){
		    m.getI().push_back(Vec3i(i,(i+1)%slice,m.getV().size()-2));
		    m.getI().push_back(Vec3i((i+1)%slice+base,i+base,m.getV().size()-1));
		}
	}
    m.applyTrans();
}
void MakeMesh::makeTorus3D(ObjMesh& m,const scalar& rad1,const scalar& rad2,const sizeType& slice1,const sizeType& slice2){
    m.getV().clear();
    for(sizeType i=0;i<slice1;i++)
    {
        scalar angI=2.0f*M_PI*(i/(scalar)slice1);
        Vec3 ctr=Vec3(cos(angI),sin(angI),0.0f)*rad1;
        Vec3 axis1=ctr*rad2/rad1;
        Vec3 axis2=Vec3::Unit(2)*rad2;
        for(sizeType j=0;j<slice2;j++)
        {
            scalar angJ=2.0f*M_PI*(j/(scalar)slice2);
            Vec3 pt=ctr+axis1*cos(angJ)+axis2*sin(angJ);
            m.getV().push_back(pt);
        }
    }
    m.getI().clear();
    for(sizeType i=0;i<slice1;i++)
    for(sizeType j=0;j<slice2;j++)
    {
        m.getI().push_back(Vec3i(GI(i  ,j  ,slice1,slice2),
                                 GI(i+1,j+1,slice1,slice2),
                                 GI(i  ,j+1,slice1,slice2)));
        m.getI().push_back(Vec3i(GI(i  ,j  ,slice1,slice2),
                                 GI(i+1,j  ,slice1,slice2),
                                 GI(i+1,j+1,slice1,slice2)));
    }
    m.applyTrans();
}
void MakeMesh::makeRing3D(ObjMesh& m,const scalar& rad1,const scalar& rad2,const scalar& rad3,const sizeType& slice){
    m.getV().clear();
    for(sizeType i=0;i<slice;i++)
    {
        scalar angI=2.0f*M_PI*(i/(scalar)slice);
        Vec3 ctr=Vec3(cos(angI),sin(angI),0.0f)*rad1;
        Vec3 axis1=ctr*rad2/rad1;
        Vec3 axis2=Vec3::Unit(2)*rad3;
        
        m.getV().push_back(ctr-axis1-axis2);
        m.getV().push_back(ctr+axis1-axis2);
        m.getV().push_back(ctr+axis1+axis2);
        m.getV().push_back(ctr-axis1+axis2);
    }
    m.getI().clear();
    for(sizeType i=0;i<slice;i++)
    for(sizeType j=0;j<4;j++)
    {
        m.getI().push_back(Vec3i(GI(i  ,j  ,slice,4),
                                 GI(i+1,j+1,slice,4),
                                 GI(i  ,j+1,slice,4)));
        m.getI().push_back(Vec3i(GI(i  ,j  ,slice,4),
                                 GI(i+1,j  ,slice,4),
                                 GI(i+1,j+1,slice,4)));
    }
    m.applyTrans();
}
void MakeMesh::makeGrid(ObjMesh& m,const Vec2i& slice){
	m.getV().clear();
	m.getI().clear();
	for(sizeType x=0;x<=slice.x();x++)
	for(sizeType y=0;y<=slice.y();y++){
		m.getV().push_back(Vec3(x/(scalar)slice.x(),y/(scalar)slice.y(),0.0f));
		if(x<slice.x() && y<slice.y()){
			m.getI().push_back(
		   Vec3i(GI(x  ,y  ,slice.x()+1,slice.y()+1),
				 GI(x+1,y  ,slice.x()+1,slice.y()+1),
				 GI(x+1,y+1,slice.x()+1,slice.y()+1)));
			m.getI().push_back(
		   Vec3i(GI(x  ,y  ,slice.x()+1,slice.y()+1),
				 GI(x+1,y+1,slice.x()+1,slice.y()+1),
				 GI(x  ,y+1,slice.x()+1,slice.y()+1)));
		}
	}
}