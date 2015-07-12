#ifndef COMPUTE_WEIGHT_INTERFACE_H
#define COMPUTE_WEIGHT_INTERFACE_H

#include "GridBasic.h"

PRJ_BEGIN

class ComputeWeightInterface
{
public:
    //interface
    virtual ~ComputeWeightInterface() {}
    virtual void computeDualCellWeight(const ScalarField& nodalSolidPhi,MACVelocityField& weight) const;
    virtual void computeCellWeight(const ScalarField& nodalSolidPhi,ScalarField& cellWeight) const;
    virtual void computeDualCellWeight3DOld(const ScalarField& nodalSolidPhi,ScalarField& weight,const sizeType& axis) const;
    virtual void computeDualCellWeight2DOld(const ScalarField& nodalSolidPhi,ScalarField& weight,const sizeType& axis) const;
    virtual void computeDualCellWeight3D(const ScalarField& nodalSolidPhi,ScalarField& weight,const sizeType& axis) const;
    virtual void computeDualCellWeight2D(const ScalarField& nodalSolidPhi,ScalarField& weight,const sizeType& axis) const;
    //helper
    static void cycle(scalar* arr,sizeType size);
    static scalar fractionCell1D(const scalar& l,const scalar& r);
    static scalar fractionCell2D(const scalar& bl,const scalar& br,const scalar& tl,const scalar& tr);
    static scalar fractionCell3D(const scalar& v0,const scalar& v1,const scalar& v2,const scalar& v3,
                                 const scalar& v4,const scalar& v5,const scalar& v6,const scalar& v7);
    static scalar fractionTet(const scalar& v0,const scalar& v1,const scalar& v2,const scalar& v3,unsigned char tag);
    static scalar fractionTet(const scalar& v0,const scalar& v1,const scalar& v2,const scalar& v3);
    static int init();
    //debug
    static scalar volTet(const Vec3& p0,const Vec3& p1,const Vec3& p2,const Vec3& p3);
    static scalar fractionTetBF(const scalar& v0,const scalar& v1,const scalar& v2,const scalar& v3,
                                const Vec3& p0,const Vec3& p1,const Vec3& p2,const Vec3& p3);
    static scalar fractionCell3DBF(const scalar& v0,const scalar& v1,const scalar& v2,const scalar& v3,
                                   const scalar& v4,const scalar& v5,const scalar& v6,const scalar& v7,
                                   const Vec3& p0,const Vec3& p1,const Vec3& p2,const Vec3& p3,
                                   const Vec3& p4,const Vec3& p5,const Vec3& p6,const Vec3& p7);
    static scalar volTri(const Vec3& p0,const Vec3& p1,const Vec3& p2);
    static scalar fractionCell2DBF(const scalar& bl,const scalar& br,const scalar& tl,const scalar& tr);
    static void debugFractionTet();
    static void debugFractionCell2D();
    static void debugFractionCell3D();
    //weight in cell
    static scalar W1,W2,W3,W4,W5,W6;
};

PRJ_END

#endif