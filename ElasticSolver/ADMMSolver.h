#ifndef ADMM_H
#define ADMM_H

#include "solvers/LinearSolver.h"

PRJ_BEGIN

template <typename KERNEL_TYPE=Kernel<scalarD> >
class MappedMat
{
public:
    typedef typename KERNEL_TYPE::Scalar T;
    typedef typename KERNEL_TYPE::Vec Vec;
    typedef typename KERNEL_TYPE::Mat Mat;
    MappedMat(const Mat& A,const vector<sizeType>& cols)
        :_A(A),_cols(cols) {}
    void mul(const Vec& x,Vec& b) const {
        b.setZero(rows());
        sizeType nrC=cols();
        for(sizeType i=0; i<nrC; i++)
            b+=_A.col(_cols[i])*x[i];
    }
    void mulT(const Vec& x,Vec& b) const {
        b.setZero(cols());
        sizeType nrC=cols();
        OMP_PARALLEL_FOR_
        for(sizeType i=0; i<nrC; i++)
            b[i]=_A.col(_cols[i]).dot(x);
    }
    void calcATA(Mat& ATA) const {
        sizeType nrC=cols();
        ATA.resize(nrC,nrC);
        OMP_PARALLEL_FOR_
        for(sizeType r=0; r<nrC; r++)
            for(sizeType c=r; c<nrC; c++)
                ATA(r,c)=ATA(c,r)=
                             _A.col(_cols[r]).dot(_A.col(_cols[c]));
    }
    sizeType rows() const {
        return _A.rows();
    }
    sizeType cols() const {
        return (sizeType)_cols.size();
    }
    T colNorm(sizeType i) const {
        return _A.col(_cols[i]).norm();
    }
protected:
    const Mat& _A;
    const vector<sizeType>& _cols;
};
template <typename KERNEL_TYPE=Kernel<scalarD> >
class ADMMKrylovMatrix : public KrylovMatrix<scalarD,KERNEL_TYPE>
{
public:
    typedef typename KERNEL_TYPE::Scalar T;
    typedef typename KERNEL_TYPE::Vec Vec;
    typedef typename KERNEL_TYPE::Mat Mat;
    ADMMKrylovMatrix(const MappedMat<KERNEL_TYPE>& A,T betaX,T betaZ)
        :_A(A),_betaX(betaX),_betaZ(betaZ),_tmp(A.rows()) {}
    virtual void multiply(const Vec& b,Vec& out) const {
        Vec& tmp=const_cast<Vec&>(_tmp);
        _A.mul(b,tmp);
        _A.mulT(tmp,out);
        KERNEL_TYPE::scale(_betaZ,out);
        KERNEL_TYPE::addScaled(_betaX,b,out);
    }
    virtual sizeType n() const {
        return _A.cols();
    }
    Vec _tmp;
    const MappedMat<KERNEL_TYPE>& _A;
    const T _betaX,_betaZ;
};
template <typename KERNEL_TYPE=Kernel<scalarD> >
class ADMMSolver
{
    typedef typename KERNEL_TYPE::Scalar T;
    typedef typename KERNEL_TYPE::Vec Vec;
    typedef typename KERNEL_TYPE::Mat Mat;
    typedef PCGSolver<T,KERNEL_TYPE,DiagonalPreconSolver<T,KERNEL_TYPE> > SOLVER;
public:
    ADMMSolver(const Matd& A):_A(A,_map) {
        setSolverParameter();
        setReport(true);
    }
    void setReport(bool report) {
        _report=report;
    }
    void setSolverParameter(T relTol=1E-3f,T relTolBound=0.1f,sizeType maxIter=1000) {
        _relTol=relTol;
        _relTolBound=relTolBound;
        _maxIter=maxIter;
        _gamma=1.618f;
    }
    void setProb(const Vec& x,const Vec& b,T p,T thres) {
        //setup Acols
        _map.clear();
        for(sizeType i=0; i<(sizeType)x.size(); i++)
            if(x[i] > thres)_map.push_back(i);
        sizeType nrC=(sizeType)_map.size();
        _x.resize(nrC);
        _w.resize(nrC);
        _b=&b;
        for(sizeType i=0; i<nrC; i++) {
            _x[i]=x[_map[i]];
            _w[i]=pow(_x[i],(T)(p-1.0f));
        }

        //setup beta
        T meanB=KERNEL_TYPE::norm1(b)/(T)b.size();
        _betaX=0.3f/meanB;
        _betaZ=3.0f/meanB;
    }
    void checkBound(T sqrDelta) const {
        Vec tmpX(_A.cols());
        KERNEL_TYPE::copy(_x,tmpX);
        KERNEL_TYPE::cwiseClamp(-numeric_limits<T>::max(),0.0f,tmpX);
        INFOV("MaxVioPos: %f",KERNEL_TYPE::absMax(tmpX))

        Vec tmpZ(_A.rows());
        _A.mul(_x,tmpZ);
        KERNEL_TYPE::sub(tmpZ,*_b,tmpZ);
        INFOV("VioBound: %f",KERNEL_TYPE::norm(tmpZ)/KERNEL_TYPE::norm(*_b))
    }
    bool solve(T sqrDelta) {
        //initialize x,r,z,w
        sizeType nrX=_A.cols();
        sizeType nrZ=_b->size();
        prepareInnerLoop();

        _Xr.resize(nrX);
        _lamX.resize(nrX);
        _tmpX.resize(nrX);
        KERNEL_TYPE::zero(_Xr);
        KERNEL_TYPE::zero(_lamX);

        _Zr.resize(nrZ);
        _lamZ.resize(nrZ);
        _tmpZ.resize(nrZ);
        KERNEL_TYPE::zero(_Zr);
        KERNEL_TYPE::zero(_lamZ);

        //the main loop
        KERNEL_TYPE::copy(_x,_xLast);
        T absTol=_relTol*_relTol*KERNEL_TYPE::dot(*_b,*_b);
        T absTolNeg=_relTolBound*_relTolBound;
        T absTolBound=_relTolBound*_relTolBound*sqrDelta;
        for(sizeType i=0; i<_maxIter; i++) {
            //update Zr
            _A.mul(_x,_tmpZ);
            KERNEL_TYPE::sub(_tmpZ,*_b,_tmpZ);
            KERNEL_TYPE::ncopy(_tmpZ,_Zr);
            KERNEL_TYPE::addScaled(1.0f/_betaZ,_lamZ,_Zr);
            T nz=KERNEL_TYPE::dot(_Zr,_Zr);
            if(nz > sqrDelta)
                KERNEL_TYPE::scale(sqrt(sqrDelta/nz),_Zr);

            //update Xr
            KERNEL_TYPE::copy(_x,_Xr);
            KERNEL_TYPE::addScaled(-1.0f/_betaX,_lamX,_Xr);
            KERNEL_TYPE::cwiseClamp(0.0f,numeric_limits<T>::max(),_Xr);

            //update x
            KERNEL_TYPE::copy(_x,_xLast);
            updateX();

            //update lambdaZ
            _A.mul(_x,_tmpZ);
            KERNEL_TYPE::sub(_tmpZ,*_b,_tmpZ);
            KERNEL_TYPE::add(_tmpZ,_Zr,_tmpZ);
            KERNEL_TYPE::addScaled(-_gamma*_betaZ,_tmpZ,_lamZ);

            //update lambdaX
            KERNEL_TYPE::sub(_x,_Xr,_tmpX);
            KERNEL_TYPE::addScaled(-_gamma*_betaX,_tmpX,_lamX);

            //stop condition
            T currNeg=KERNEL_TYPE::dot(_tmpX,_tmpX);
            T currBound=KERNEL_TYPE::dot(_tmpZ,_tmpZ);
            KERNEL_TYPE::sub(_x,_xLast,_tmpX);
            T curr=KERNEL_TYPE::dot(_tmpX,_tmpX);
            if(curr < absTol && currNeg < absTolNeg && currBound < absTolBound)
                return true;
            //report
            if(_report)
                INFOV("Iter%d: RelTol:(%f,%f) VioNeg:(%f,%f) VioBound:(%f,%f)",
                      i,curr,absTol,currNeg,absTolNeg,currBound,absTolBound)
            }
        return false;
    }
    void fetchX(Vec& x) const {
        x.setZero();
        OMP_PARALLEL_FOR_
        for(sizeType i=0; i<(sizeType)_map.size(); i++)
            x[_map[i]]=_x[i];
    }
protected:
    void prepareInnerLoop() {
#define THRES_ADMM_DIRECT 500000	//so that we never use Krylov method
        if(_A.cols() > THRES_ADMM_DIRECT) {
            _sol.reset(new SOLVER());
            _sol->setSolverParameters(1E-8f,_maxIter);
            _sol->setKrylovMatrix(boost::shared_ptr<KrylovMatrix<T,KERNEL_TYPE> >(new ADMMKrylovMatrix<KERNEL_TYPE>(_A,_betaX,_betaZ)));
            DiagonalPreconSolver<T,KERNEL_TYPE>* pre=((DiagonalPreconSolver<T,KERNEL_TYPE>*)_sol->getPre());

            sizeType nrX=_A.cols();
            Vec& diag=pre->_diag;
            diag.resize(nrX);
            for(sizeType i=0; i<nrX; i++) {
                T d=_A.colNorm(i);
                diag[i]=1.0f/(_betaX+_betaZ*d*d);
            }
        } else {
            _sol.reset((SOLVER*)NULL);
            Mat invA;
            _A.calcATA(invA);
            invA*=_betaZ;
            invA.diagonal().array()+=_betaX;
            _solD=invA.llt();
        }
#undef THRES_ADMM_DIRECT
    }
    void updateX() {
        KERNEL_TYPE::sub(*_b,_Zr,_tmpZ);
        KERNEL_TYPE::scale(_betaZ,_tmpZ);
        KERNEL_TYPE::add(_lamZ,_tmpZ,_tmpZ);
        _A.mulT(_tmpZ,_tmpX);
        KERNEL_TYPE::addScaled(_betaX,_Xr,_tmpX);
        KERNEL_TYPE::add(_lamX,_tmpX,_tmpX);
        KERNEL_TYPE::sub(_tmpX,_w,_tmpX);
        if(_sol)_sol->solve(_tmpX,_x);
        else _x=_solD.solve(_tmpX);
    }
private:
    //data
    const MappedMat<KERNEL_TYPE> _A;
    Vec _x,_w;
    const Vec *_b;
    vector<sizeType> _map;
    //krylov solver
    boost::shared_ptr<SOLVER> _sol;
    Eigen::LLT<Mat> _solD;
    //stopping param
    T _relTol,_relTolBound,_gamma;
    sizeType _maxIter;
    //runtime param
    T _betaZ,_betaX;
    Vec _tmpZ,_tmpX;
    Vec _lamZ,_lamX;
    Vec _Zr,_Xr;
    Vec _xLast;
    bool _report;
};

PRJ_END

#endif