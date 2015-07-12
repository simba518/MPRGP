#ifndef FEM_RS_REDUCED_SYSTEM_H
#define FEM_RS_REDUCED_SYSTEM_H

#include "FEMRigidReducedSystem.h"

PRJ_BEGIN

struct PolyCoefEfficient;
class CoupledRSReducedMCalculator : public ReducedMCalculator
{
public:
    CoupledRSReducedMCalculator(FEMReducedSystem& sys);
    virtual void precompute();
    virtual void init();
    virtual void calculateM(const Vec& S,Matd& M,const Vec& MRHSL,Vec& RHS,sizeType r,sizeType c,FEMSystemMatrix* U);
    virtual Vec calculateFE(const Vec& S,Matd& K);
    //utility
    virtual void getPosL(Vec& L) const;
    virtual Vec LtoF(const Vec& L) const;
    virtual void setPos(const Vec& P);
    virtual sizeType sizeL() const;
protected:
    virtual void setRSConstraint(Vec& X) const;
    virtual void calculateMX(const Vec& P);
    void buildMatrix(Eigen::SparseMatrix<scalarD,0,sizeType>& PRJ,bool cons);
    void debugDLDS();
    //data
    const boost::unordered_map<sizeType,Vec3>& _constraints;
    Eigen::SparseMatrix<scalarD,0,sizeType> _G,_V;
    Eigen::SparseLU<Eigen::SparseMatrix<scalarD,0,sizeType> > _invGTVG;
    //temporary data
    Vec _lastS,_L;
    Matd _GU,_DLDS;
};
class CRCoupledRSReducedMCalculatorBase : public CoupledRSReducedMCalculator
{
public:
    CRCoupledRSReducedMCalculatorBase(FEMReducedSystem& sys);
    virtual void buildBasis();
protected:
    Matd calculateFtoB();
    void debugRSB1(const Vec& col,sizeType i);
    void debugRSB2(const Vec& col,sizeType i,sizeType j);
};
class CRCoupledRSReducedMCalculatorRef : public CRCoupledRSReducedMCalculatorBase
{
public:
    CRCoupledRSReducedMCalculatorRef(FEMReducedSystem& sys);
    virtual void precompute();
    virtual void calculateM(const Vec& S,Matd& M,const Vec& MRHSL,Vec& RHS,sizeType r,sizeType c,FEMSystemMatrix* U);
    virtual Vec calculateFE(const Vec& S,Matd& K);
    virtual Vec LtoF(const Vec& L) const;
    virtual sizeType sizeL() const;
protected:
    virtual void calculateMX(const Vec& P);
    virtual const Matd& getBasis() const;
    Matd _matM,_BCub;
    boost::shared_ptr<Matd> _BTMB;
};
class CRCoupledRSReducedMCalculator : public CRCoupledRSReducedMCalculatorRef
{
public:
    CRCoupledRSReducedMCalculator(FEMReducedSystem& sys);
    virtual void precompute();
protected:
    virtual const Matd& getBasis() const;
    Matd _BInvLT;
};
class TRCoupledRSReducedMCalculator : public CRCoupledRSReducedMCalculator
{
public:
    TRCoupledRSReducedMCalculator(FEMReducedSystem& sys);
    virtual void precompute();
protected:
    virtual void calculateMX(const Vec& P);
    boost::shared_ptr<PolyCoefEfficient> _coefP,_coefQ;
};

PRJ_END

#endif