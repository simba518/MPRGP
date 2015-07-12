#ifndef FEM_WARPPED_REDUCED_SYSTEM_H
#define FEM_WARPPED_REDUCED_SYSTEM_H

#include "FEMReducedSystem.h"

PRJ_BEGIN

class WrappedReducedMCalculator : public ReducedMCalculator
{
public:
    WrappedReducedMCalculator(FEMReducedSystem& sys);
    virtual void precompute();
    virtual void calculateM(const Vec& S,Matd& M,const Vec& MRHSL,Vec& RHS,sizeType r,sizeType c,FEMSystemMatrix* U);
    virtual Vec calculateFE(const Vec& S,Matd& K);
    //utility
    virtual void setPos(const Vec& P);
private:
    Eigen::SparseMatrix<scalarD,0,sizeType> _R;
    Matd _Gamma,_UDense;
};

PRJ_END

#endif