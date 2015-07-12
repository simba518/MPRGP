#include "ComputeWeightInterface.h"

USE_PRJ_NAMESPACE

void ComputeWeightInterface::computeDualCellWeight(const ScalarField& nodalSolidPhi,MACVelocityField& weight) const
{
    /*
    //old code
    if(nodalSolidPhi.getDim() == 3)
    {
    	computeDualCellWeight3DOld(nodalSolidPhi,weight.getGu(),0);
    	computeDualCellWeight3DOld(nodalSolidPhi,weight.getGv(),1);
    	computeDualCellWeight3DOld(nodalSolidPhi,weight.getGw(),2);
    }
    else if(nodalSolidPhi.getDim() == 2)
    {
    	computeDualCellWeight2DOld(nodalSolidPhi,weight.getGu(),0);
    	computeDualCellWeight2DOld(nodalSolidPhi,weight.getGv(),1);
    }
    else ASSERT(false);
    */

    if(nodalSolidPhi.getDim() == 3) {
        computeDualCellWeight3D(nodalSolidPhi,weight.getGu(),0);
        computeDualCellWeight3D(nodalSolidPhi,weight.getGv(),1);
        computeDualCellWeight3D(nodalSolidPhi,weight.getGw(),2);
    } else if(nodalSolidPhi.getDim() == 2) {
        computeDualCellWeight2D(nodalSolidPhi,weight.getGu(),0);
        computeDualCellWeight2D(nodalSolidPhi,weight.getGv(),1);
    } else ASSERT(false)
    }
void ComputeWeightInterface::computeCellWeight(const ScalarField& nodalSolidPhi,ScalarField& cellWeight) const
{
    ASSERT(nodalSolidPhi.getNrCell() == cellWeight.getNrCell())
    ASSERT(!nodalSolidPhi.isCenter() && cellWeight.isCenter())
    if(nodalSolidPhi.getDim() == 3) {
        const Vec3i nrCell=cellWeight.getNrCell();
        for(sizeType x=0; x<nrCell.x(); x++)
            for(sizeType y=0; y<nrCell.y(); y++)
                for(sizeType z=0; z<nrCell.z(); z++) {
                    cellWeight.get(Vec3i(x,y,z))=
                        1.0f-fractionCell3D(nodalSolidPhi.get(Vec3i(x  ,y  ,z  )),
                                            nodalSolidPhi.get(Vec3i(x  ,y  ,z+1)),
                                            nodalSolidPhi.get(Vec3i(x+1,y  ,z+1)),
                                            nodalSolidPhi.get(Vec3i(x+1,y  ,z  )),

                                            nodalSolidPhi.get(Vec3i(x  ,y+1,z  )),
                                            nodalSolidPhi.get(Vec3i(x  ,y+1,z+1)),
                                            nodalSolidPhi.get(Vec3i(x+1,y+1,z+1)),
                                            nodalSolidPhi.get(Vec3i(x+1,y+1,z  )));
                }
    } else if(nodalSolidPhi.getDim() == 2) {
        const Vec3i nrCell=cellWeight.getNrCell();
        for(sizeType x=0; x<nrCell.x(); x++)
            for(sizeType y=0; y<nrCell.y(); y++) {
                cellWeight.get(Vec3i(x,y,0))=
                    1.0f-fractionCell2D(nodalSolidPhi.get(Vec3i(x  ,y  ,0)),
                                        nodalSolidPhi.get(Vec3i(x+1,y  ,0)),
                                        nodalSolidPhi.get(Vec3i(x  ,y+1,0)),
                                        nodalSolidPhi.get(Vec3i(x+1,y+1,0)));
            }
    } else ASSERT(false);
}
void ComputeWeightInterface::computeDualCellWeight3DOld(const ScalarField& N,ScalarField& weight,const sizeType& axis) const
{
    Vec3i a0,a1,a2;
    if(axis == 0) {
        a0=Vec3i(0,1,0);
        a1=Vec3i(0,0,1);
        a2=Vec3i(0,1,1);
    } else if(axis == 1) {
        a0=Vec3i(1,0,0);
        a1=Vec3i(0,0,1);
        a2=Vec3i(1,0,1);
    } else {
        a0=Vec3i(1,0,0);
        a1=Vec3i(0,1,0);
        a2=Vec3i(1,1,0);
    }

    const Vec3i nrPoint=weight.getNrPoint();
    OMP_PARALLEL_FOR_
    for(sizeType x=0; x<nrPoint.x(); x++)
        for(sizeType y=0; y<nrPoint.y(); y++)
            for(sizeType z=0; z<nrPoint.z(); z++) {
                const Vec3i id(x,y,z);
                weight.get(id)=1.0f-fractionCell2D(N.get(id),N.get(id+a0),N.get(id+a1),N.get(id+a2));
                //weight.get(id)=min<scalar>(max<scalar>(weight.get(id),0.0f),1.0f);
            }
}
void ComputeWeightInterface::computeDualCellWeight2DOld(const ScalarField& N,ScalarField& weight,const sizeType& axis) const
{
    Vec3i a0;
    if(axis == 0)
        a0=Vec3i(0,1,0);
    else
        a0=Vec3i(1,0,0);

    const Vec3i nrPoint=weight.getNrPoint();
    OMP_PARALLEL_FOR_
    for(sizeType x=0; x<nrPoint.x(); x++)
        for(sizeType y=0; y<nrPoint.y(); y++) {
            const Vec3i id(x,y,0);
            weight.get(id)=1.0f-fractionCell1D(N.getSafe(id),N.getSafe(id+a0));
            //weight.get(id)=min<scalar>(max<scalar>(weight.get(id),0.0f),1.0f);
        }
}
void ComputeWeightInterface::computeDualCellWeight3D(const ScalarField& N,ScalarField& weight,const sizeType& axis) const
{
    const Vec3i unit=Vec3i::Unit(axis);
    const Vec3i nrPoint=weight.getNrPoint();
    weight.init(0.0f);
    for(sizeType x=0; x<nrPoint.x(); x++)
        for(sizeType y=0; y<nrPoint.y(); y++)
            for(sizeType z=0; z<nrPoint.z(); z++) {
                const Vec3i A=Vec3i(x,y,z),B=A-unit;
                weight.get(A)=
                    1.0f-fractionCell3D((N.getSafe(A+Vec3i(0,0,0))+N.getSafe(B+Vec3i(0,0,0)))/2.0f,
                                        (N.getSafe(A+Vec3i(0,0,1))+N.getSafe(B+Vec3i(0,0,1)))/2.0f,
                                        (N.getSafe(A+Vec3i(1,0,1))+N.getSafe(B+Vec3i(1,0,1)))/2.0f,
                                        (N.getSafe(A+Vec3i(1,0,0))+N.getSafe(B+Vec3i(1,0,0)))/2.0f,

                                        (N.getSafe(A+Vec3i(0,1,0))+N.getSafe(B+Vec3i(0,1,0)))/2.0f,
                                        (N.getSafe(A+Vec3i(0,1,1))+N.getSafe(B+Vec3i(0,1,1)))/2.0f,
                                        (N.getSafe(A+Vec3i(1,1,1))+N.getSafe(B+Vec3i(1,1,1)))/2.0f,
                                        (N.getSafe(A+Vec3i(1,1,0))+N.getSafe(B+Vec3i(1,1,0)))/2.0f);
            }
}
void ComputeWeightInterface::computeDualCellWeight2D(const ScalarField& N,ScalarField& weight,const sizeType& axis) const
{
    const Vec3i unit=Vec3i::Unit(axis);
    const Vec3i nrPoint=weight.getNrPoint();
    weight.init(0.0f);
    for(sizeType x=0; x<nrPoint.x(); x++)
        for(sizeType y=0; y<nrPoint.y(); y++) {
            Vec3i A(x,y,0),B=A-unit;
            weight.get(A)=
                1.0f-fractionCell2D((N.getSafe(A+Vec3i(0,0,0))+N.getSafe(B+Vec3i(0,0,0)))/2.0f,
                                    (N.getSafe(A+Vec3i(1,0,0))+N.getSafe(B+Vec3i(1,0,0)))/2.0f,
                                    (N.getSafe(A+Vec3i(0,1,0))+N.getSafe(B+Vec3i(0,1,0)))/2.0f,
                                    (N.getSafe(A+Vec3i(1,1,0))+N.getSafe(B+Vec3i(1,1,0)))/2.0f);
        }
}

void ComputeWeightInterface::cycle(scalar* arr,sizeType size)
{
    scalar t=arr[0];
    for(sizeType i=0; i<size-1; ++i)
        arr[i]=arr[i+1];
    arr[size-1]=t;
}
scalar ComputeWeightInterface::fractionCell1D(const scalar& l,const scalar& r)
{
    if(l < 0.0f && r < 0.0f)
        return 1.0f;
    if(l < 0.0f && r >= 0.0f)
        return l/(l-r);
    if(l >= 0.0f && r < 0.0f)
        return r/(r-l);
    else
        return 0.0f;
}
scalar ComputeWeightInterface::fractionCell2D(const scalar& bl,const scalar& br,const scalar& tl,const scalar& tr)
{
    sizeType insideCount=(bl<0.0f?1:0)+(tl<0.0f?1:0)+(br<0.0f?1:0)+(tr<0.0f?1:0);
    scalar list[]= {bl,br,tr,tl};

    if(insideCount == 4) {
        return 1.0f;
    } else if(insideCount == 3) {
        //rotate until the positive value is in the first position
        while(list[0]<0.0f)
            cycle(list,4);

        //Work out the area of the exterior triangle
        scalar side0=1.0f-fractionCell1D(list[0],list[3]);
        scalar side1=1.0f-fractionCell1D(list[0],list[1]);
        return 1.0f-0.5f*side0*side1;
    } else if(insideCount == 2) {
        //rotate until a negative value is in the first position, and the next negative is in either slot 1 or 2.
        while(list[0] >= 0.0f || !(list[1] < 0.0f || list[2] < 0.0f))
            cycle(list,4);

        if(list[1] < 0.0f) {
            //the matching signs are adjacent
            scalar sideLeft=fractionCell1D(list[0],list[3]);
            scalar sideRight=fractionCell1D(list[1],list[2]);
            return 0.5f*(sideLeft+sideRight);
        } else {
            //matching signs are diagonally opposite
            //determine the centre point's sign to disambiguate this case
            scalar middlePoint=0.25f*(list[0]+list[1]+list[2]+list[3]);
            if(middlePoint < 0.0f) {
                scalar area=0.0f;

                //first triangle (top left)
                scalar side1=1.0f-fractionCell1D(list[0],list[3]);
                scalar side3=1.0f-fractionCell1D(list[2],list[3]);

                area+=0.5f*side1*side3;

                //second triangle (top right)
                scalar side2=1.0f-fractionCell1D(list[2],list[1]);
                scalar side0=1.0f-fractionCell1D(list[0],list[1]);
                area+=0.5f*side0*side2;

                return 1.0f-area;
            } else {
                scalar area=0.0f;

                //first triangle (bottom left)
                scalar side0=fractionCell1D(list[0],list[1]);
                scalar side1=fractionCell1D(list[0],list[3]);
                area+=0.5f*side0*side1;

                //second triangle (top right)
                scalar side2=fractionCell1D(list[2],list[1]);
                scalar side3=fractionCell1D(list[2],list[3]);
                area+=0.5f*side2*side3;
                return area;
            }
        }
    } else if(insideCount == 1) {
        //rotate until the negative value is in the first position
        while(list[0] >= 0.0f)
            cycle(list,4);

        //Work out the area of the interior triangle, and subtract from 1.
        scalar side0=fractionCell1D(list[0],list[3]);
        scalar side1=fractionCell1D(list[0],list[1]);
        return 0.5f*side0*side1;
    } else {
        return 0.0f;
    }
}
scalar ComputeWeightInterface::fractionCell3D(const scalar& v0,const scalar& v1,const scalar& v2,const scalar& v3,
        const scalar& v4,const scalar& v5,const scalar& v6,const scalar& v7)
{
    static int verbose=init();
    return fractionTet(v2,v0,v4,v3)*W1+
           fractionTet(v2,v7,v3,v4)*W2+
           fractionTet(v2,v4,v6,v7)*W3+
           fractionTet(v2,v4,v5,v6)*W4+
           fractionTet(v2,v1,v5,v4)*W5+
           fractionTet(v2,v0,v1,v4)*W6;
}

scalar ComputeWeightInterface::fractionTet(const scalar& v0,const scalar& v1,const scalar& v2,const scalar& v3,unsigned char tag)
{
    switch(tag) {
    case 0:
        return 0.0f;
    case 1:
        return -v0*v0*v0/std::max((v1-v0)*(v2-v0)*(v3-v0),EPS);
    case 2:
        return -v1*v1*v1/std::max((v0-v1)*(v2-v1)*(v3-v1),EPS);
    case 3:
        return v1*v1/std::max((v2-v1)*(v3-v1),EPS)+
               v0*v2*v1/std::max((v2-v0)*(v2-v1)*(v3-v1),EPS)+
               v0*v0*v3/std::max((v3-v0)*(v2-v0)*(v3-v1),EPS);
    case 4:
        return -v2*v2*v2/std::max((v0-v2)*(v1-v2)*(v3-v2),EPS);
    case 5:
        return v2*v2/std::max((v3-v2)*(v1-v2),EPS)+
               v0*v3*v2/std::max((v3-v0)*(v3-v2)*(v1-v2),EPS)+
               v0*v0*v1/std::max((v1-v0)*(v3-v0)*(v1-v2),EPS);
    case 6:
        return v1*v1/std::max((v3-v1)*(v0-v1),EPS)+
               v2*v3*v1/std::max((v3-v2)*(v3-v1)*(v0-v1),EPS)+
               v2*v2*v0/std::max((v0-v2)*(v3-v2)*(v0-v1),EPS);
    case 8:
        return -v3*v3*v3/std::max((v0-v3)*(v1-v3)*(v2-v3),EPS);
    case 15:
        return 1.0f;
    }
    return 1.0f-fractionTet(-v0,-v1,-v2,-v3,15-tag);
}
scalar ComputeWeightInterface::fractionTet(const scalar& v0,const scalar& v1,const scalar& v2,const scalar& v3)
{
    unsigned char tag=(v0 < 0.0f ? 1 : 0)+
                      (v1 < 0.0f ? 2 : 0)+
                      (v2 < 0.0f ? 4 : 0)+
                      (v3 < 0.0f ? 8 : 0);
    return fractionTet(v0,v1,v2,v3,tag);
}
int ComputeWeightInterface::init()
{
    Vec3 v0(0.0f,0.0f,0.0f);
    Vec3 v1(1.0f,0.0f,0.0f);
    Vec3 v2(1.0f,0.0f,1.0f);
    Vec3 v3(0.0f,0.0f,1.0f);

    Vec3 v4(0.0f,1.0f,0.0f);
    Vec3 v5(1.0f,1.0f,0.0f);
    Vec3 v6(1.0f,1.0f,1.0f);
    Vec3 v7(0.0f,1.0f,1.0f);

    W1=volTet(v2,v0,v4,v3);
    W2=volTet(v2,v7,v3,v4);
    W3=volTet(v2,v4,v6,v7);
    W4=volTet(v2,v4,v5,v6);
    W5=volTet(v2,v1,v5,v4);
    W6=volTet(v2,v0,v1,v4);

    scalar TOTAL=W1+W2+W3+W4+W5+W6;
    W1/=TOTAL;
    W2/=TOTAL;
    W3/=TOTAL;
    W4/=TOTAL;
    W5/=TOTAL;
    W6/=TOTAL;
    return 0;
}

scalar ComputeWeightInterface::volTet(const Vec3& p0,const Vec3& p1,const Vec3& p2,const Vec3& p3)
{
    return std::abs((p0-p3).dot((p1-p3).cross(p2-p3)))/6.0f;
}
scalar ComputeWeightInterface::fractionTetBF(const scalar& v0,const scalar& v1,const scalar& v2,const scalar& v3,
        const Vec3& p0,const Vec3& p1,const Vec3& p2,const Vec3& p3)
{
#define interp(va,vb,pa,pb) ((pa*vb-pb*va)/(vb-va))
    const scalar coef=1.0f/volTet(p0,p1,p2,p3);
    unsigned char tag=(v0 < 0.0f ? 1 : 0)+
                      (v1 < 0.0f ? 2 : 0)+
                      (v2 < 0.0f ? 4 : 0)+
                      (v3 < 0.0f ? 8 : 0);
    switch(tag) {
    case 0:
        return 0.0f;
    case 1:
        return coef*(volTet(p0,interp(v0,v1,p0,p1),interp(v0,v2,p0,p2),interp(v0,v3,p0,p3)));
    case 2:
        return coef*(volTet(p1,interp(v1,v0,p1,p0),interp(v1,v2,p1,p2),interp(v1,v3,p1,p3)));
    case 3:
        return coef*(volTet(p1,p0,interp(v1,v2,p1,p2),interp(v1,v3,p1,p3))+
                     volTet(p0,interp(v0,v2,p0,p2),interp(v1,v3,p1,p3),interp(v1,v2,p1,p2))+
                     volTet(p0,interp(v0,v2,p0,p2),interp(v0,v3,p0,p3),interp(v1,v3,p1,p3)));
    case 4:
        return coef*(volTet(p2,interp(v2,v0,p2,p0),interp(v2,v1,p2,p1),interp(v2,v3,p2,p3)));
    case 5:
        return coef*(volTet(p0,p2,interp(v2,v1,p2,p1),interp(v2,v3,p2,p3))+
                     volTet(p0,interp(v2,v1,p2,p1),interp(v2,v3,p2,p3),interp(v0,v1,p0,p1))+
                     volTet(p0,interp(v2,v3,p2,p3),interp(v0,v3,p0,p3),interp(v0,v1,p0,p1)));
    case 6:
        return coef*(volTet(p2,p1,interp(v1,v0,p1,p0),interp(v1,v3,p1,p3))+
                     volTet(p2,interp(v2,v3,p2,p3),interp(v1,v3,p1,p3),interp(v1,v0,p1,p0))+
                     volTet(p2,interp(v2,v3,p2,p3),interp(v2,v0,p2,p0),interp(v1,v0,p1,p0)));
    case 8:
        return coef*(volTet(p3,interp(v3,v0,p3,p0),interp(v3,v1,p3,p1),interp(v3,v2,p3,p2)));
    case 15:
        return 1.0f;
    }
    return 1.0f-fractionTetBF(-v0,-v1,-v2,-v3,p0,p1,p2,p3);
}
scalar ComputeWeightInterface::fractionCell3DBF(const scalar& v0,const scalar& v1,const scalar& v2,const scalar& v3,
        const scalar& v4,const scalar& v5,const scalar& v6,const scalar& v7,
        const Vec3& p0,const Vec3& p1,const Vec3& p2,const Vec3& p3,
        const Vec3& p4,const Vec3& p5,const Vec3& p6,const Vec3& p7)
{
    static int verbose=init();
    return fractionTetBF(v2,v0,v4,v3, p2,p0,p4,p3)*W1+
           fractionTetBF(v2,v7,v3,v4, p2,p7,p3,p4)*W2+
           fractionTetBF(v2,v4,v6,v7, p2,p4,p6,p7)*W3+
           fractionTetBF(v2,v4,v5,v6, p2,p4,p5,p6)*W4+
           fractionTetBF(v2,v1,v5,v4, p2,p1,p5,p4)*W5+
           fractionTetBF(v2,v0,v1,v4, p2,p0,p1,p4)*W6;
}
scalar ComputeWeightInterface::volTri(const Vec3& p0,const Vec3& p1,const Vec3& p2)
{
    return (p1-p0).cross(p2-p0).norm()/2.0f;
}
scalar ComputeWeightInterface::fractionCell2DBF(const scalar& v0,const scalar& v1,const scalar& v2,const scalar& v3)
{
    const Vec3 p0(0.0f,0.0f,0.0f);
    const Vec3 p1(1.0f,0.0f,0.0f);
    const Vec3 p2(0.0f,1.0f,0.0f);
    const Vec3 p3(1.0f,1.0f,0.0f);
    unsigned char tag=(v0 < 0.0f ? 1 : 0)+
                      (v1 < 0.0f ? 2 : 0)+
                      (v2 < 0.0f ? 4 : 0)+
                      (v3 < 0.0f ? 8 : 0);
    scalar sum=v0+v1+v2+v3;
    INFOV("Sum: %f",sum)
    switch(tag) {
    case 0:
        return 0.0f;
    case 1:
        return volTri(p0,interp(v0,v1,p0,p1),interp(v0,v2,p0,p2));
    case 2:
        return volTri(p1,interp(v0,v1,p0,p1),interp(v1,v3,p1,p3));
    case 3:
        return volTri(p0,p1,interp(v0,v2,p0,p2))+volTri(p0,p1,interp(v1,v3,p1,p3));
    case 4:
        return volTri(p2,interp(v0,v2,p0,p2),interp(v2,v3,p2,p3));
    case 5:
        return volTri(p0,p2,interp(v0,v1,p0,p1))+volTri(p0,p2,interp(v2,v3,p2,p3));
    case 6:
        if(sum < 0.0f)
            return 1.0f-volTri(p0,interp(v0,v1,p0,p1),interp(v0,v2,p0,p2))-
                   volTri(p3,interp(v3,v1,p3,p1),interp(v3,v2,p3,p2));
        else
            return volTri(p1,interp(v0,v1,p0,p1),interp(v1,v3,p1,p3))+
                   volTri(p2,interp(v0,v2,p0,p2),interp(v2,v3,p2,p3));
    case 8:
        return volTri(p3,interp(v1,v3,p1,p3),interp(v2,v3,p2,p3));
    case 15:
        return 1.0f;
    }
    return 1.0f-fractionCell2DBF(-v0,-v1,-v2,-v3);
}
void ComputeWeightInterface::debugFractionTet()
{
    Vec3 a(-1.0f,-3.0f,-5.0f);
    Vec3 b(1.0f,2.0f,0.0f);
    Vec3 c(2.0f,1.0f,0.0f);
    Vec3 d(5.5f,0.0f,1.0f);

    a.setRandom();
    b.setRandom();
    c.setRandom();
    d.setRandom();

    scalar va;
    scalar vb;
    scalar vc;
    scalar vd;

    for(sizeType i=0; i<16; i++) {
        va=(i&1) ? -1.0f : 1.0f;
        vb=(i&2) ? -1.0f : 1.0f;
        vc=(i&4) ? -1.0f : 1.0f;
        vd=(i&8) ? -1.0f : 1.0f;

        va*=(scalar)rand()/(scalar)RAND_MAX;
        vb*=(scalar)rand()/(scalar)RAND_MAX;
        vc*=(scalar)rand()/(scalar)RAND_MAX;
        vd*=(scalar)rand()/(scalar)RAND_MAX;

        scalar fracFBA=fractionTetBF(va,vb,vc,vd,a,b,c,d);
        scalar fracFBB=fractionTetBF(va,vc,vb,vd,a,c,b,d);
        scalar fracFBC=fractionTetBF(vd,va,vb,vc,d,a,b,c);
        scalar fracFBD=fractionTetBF(va,vb,vd,vc,a,b,d,c);

        scalar fracA=fractionTet(va,vb,vc,vd);
        scalar fracB=fractionTet(va,vd,vb,vc);
        scalar fracC=fractionTet(va,vc,vd,vb);
        scalar fracD=fractionTet(vd,va,vc,vb);

        INFOV("Fraction Tet %d:",i)
        INFOV("\tBruteForce: %f,%f,%f,%f",fracFBA,fracFBB,fracFBC,fracFBD)
        INFOV("\tFastCalcul: %f,%f,%f,%f",fracA,fracB,fracC,fracD)

        ASSERT(std::abs(fracFBA-fracFBB) < EPS)
        ASSERT(std::abs(fracFBA-fracFBC) < EPS)
        ASSERT(std::abs(fracFBA-fracFBD) < EPS)

        ASSERT(std::abs(fracFBA-fracA) < EPS)
        ASSERT(std::abs(fracFBA-fracB) < EPS)
        ASSERT(std::abs(fracFBA-fracD) < EPS)
        ASSERT(std::abs(fracFBA-fracC) < EPS)
    }
}
void ComputeWeightInterface::debugFractionCell2D()
{
    scalar va;
    scalar vb;
    scalar vc;
    scalar vd;

    for(sizeType i=0; i<16; i++) {
        va=(i&1) ? -1.0f : 1.0f;
        vb=(i&2) ? -1.0f : 1.0f;
        vc=(i&4) ? -1.0f : 1.0f;
        vd=(i&8) ? -1.0f : 1.0f;

        va*=(scalar)rand()/(scalar)RAND_MAX;
        vb*=(scalar)rand()/(scalar)RAND_MAX;
        vc*=(scalar)rand()/(scalar)RAND_MAX;
        vd*=(scalar)rand()/(scalar)RAND_MAX;

        scalar fracA=fractionCell2D(va,vb,vc,vd);
        scalar fracB=fractionCell2DBF(va,vb,vc,vd);

        INFOV("Fraction Cell 2D: %f %f",fracA,fracB)
        ASSERT(std::abs(fracA-fracB) < EPS)
    }
}
void ComputeWeightInterface::debugFractionCell3D()
{
    scalar cx=1.0f;
    scalar cy=6.0f;
    scalar cz=3.0f;

    Vec3 v0(0.0f,0.0f,0.0f);
    Vec3 v1(cx  ,0.0f,0.0f);
    Vec3 v2(cx  ,0.0f,cz  );
    Vec3 v3(0.0f,0.0f,cz  );

    Vec3 v4(0.0f,cy  ,0.0f);
    Vec3 v5(cx  ,cy  ,0.0f);
    Vec3 v6(cx  ,cy  ,cz  );
    Vec3 v7(0.0f,cy  ,cz  );

    INFOV("FractionCell3DFastCalcul: %f",fractionCell3D  (2.0f,2.0f,-2.0f,-2.0f,6.0f,6.0f, 2.0f, 2.0f))
    INFOV("FractionCell3DBruteForce: %f",fractionCell3DBF(2.0f,2.0f,-2.0f,-2.0f,6.0f,6.0f, 2.0f, 2.0f, v0,v1,v2,v3,v4,v5,v6,v7))
}

scalar ComputeWeightInterface::W1;
scalar ComputeWeightInterface::W2;
scalar ComputeWeightInterface::W3;
scalar ComputeWeightInterface::W4;
scalar ComputeWeightInterface::W5;
scalar ComputeWeightInterface::W6;