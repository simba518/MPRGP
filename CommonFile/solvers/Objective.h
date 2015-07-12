#ifndef OBJECTIVE_H
#define OBJECTIVE_H

#include "MathBasic.h"
#include "MatVec.h"

PRJ_BEGIN

template < typename T, typename KERNEL_TYPE=Kernel<T> >
struct Objective {
public:
    typedef typename KERNEL_TYPE::Vec Vec;
	//for L-BFGS Optimization
    virtual int operator()(const Vec& x,T& FX,Vec& DFDX,const T& step,bool wantGradient) {
        ASSERT_MSG(false,"Not Implemented: Objective Function Without Hessian!")
        return -1;
    }
	//For LM Optimization
    virtual int operator()(const Vec& x,Vec& fvec) {
        ASSERT_MSG(false,"Not Implemented: Least Square Function!")
        return -1;
    }
    virtual int df(const Vec& x,Eigen::MatrixXd& fjac) {
        ASSERT_MSG(false,"Not Implemented: Least Square Jacobi!")
        return -1;
    }
	//For Augmented Lagrangian Optimization
	virtual int operator()(const Vec& x,Vec& cvec,std::vector<Eigen::Triplet<T,sizeType> >& cjac) {
		return -1;
	}
	//Dimension Info
    virtual int inputs() const {
        return 0;
    }
    virtual int values() const {
        return 0;
    }
	virtual int constraints() const {
		return 0;
	}
    virtual void profileLineSearch(const sizeType& k,const Vec& x,const Vec& d,const T& step) {
        return;
    }
};

PRJ_END

#endif