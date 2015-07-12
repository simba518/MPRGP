#ifndef FEM_RIGID_REDUCED_SYSTEM_H
#define FEM_RIGID_REDUCED_SYSTEM_H

#include "FEMReducedSystem.h"

PRJ_BEGIN

class RigidReducedMCalculator : public ReducedMCalculator
{
public:
    RigidReducedMCalculator(FEMReducedSystem& sys,boost::shared_ptr<ReducedMCalculator> innerMCalc);
    void setRigid(bool rigid);
    virtual void calculateM(const Vec& S,Matd& M,const Vec& MRHSL,Vec& RHS,sizeType r,sizeType c,FEMSystemMatrix* U);
    virtual Vec calculateFE(const Vec& S,Matd& K);
    virtual void getPos(Vec& P) const;
    virtual void getPosL(Vec& L) const;
    virtual Vec LtoF(const Vec& L) const;
    virtual void setPos(const Vec& P);
    virtual sizeType size() const;
    virtual sizeType sizeL() const;
protected:
    Matd _UDense;
    boost::shared_ptr<ReducedMCalculator> _innerMCalc;
};

PRJ_END

#endif