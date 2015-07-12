#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include "MathBasic.h"
#include "KrylovMatrix.h"
#include "Callback.h"
#include <boost/shared_ptr.hpp>

PRJ_BEGIN

//============================================================================
// Linear Solver Interface
template <typename T,typename KERNEL_TYPE>
struct Solver {
public:
    enum SOLVER_RESULT {
        SUCCESSFUL,
        NOT_CONVERGENT,
        USER_REQUEST,
        NEGATIVE_EIGENVALUE,
    };
public:
    typedef typename KERNEL_TYPE::Vec Vec;
    Solver():_cb((Callback<T,KERNEL_TYPE>*)NULL) {}
    virtual ~Solver() {}
    virtual void setMatrix(const FixedSparseMatrix<T,KERNEL_TYPE>& matrix,bool syncPrecon) =0;
    virtual void setKrylovMatrix(boost::shared_ptr<KrylovMatrix<T,KERNEL_TYPE> > matrix){ASSERT_MSG(false,"Not Implemented!");}
    virtual void setSolverParameters(T toleranceFactor,sizeType maxIterations){}
    virtual void resize(const sizeType& n){}
    virtual Solver<T,KERNEL_TYPE>* getPre() {return NULL;}
    virtual const Solver<T,KERNEL_TYPE>* getPre() const {return NULL;}
    virtual SOLVER_RESULT solve(const Vec& rhs,Vec& result) =0;
    SOLVER_RESULT solveMatrix(const SparseMatrix<T,KERNEL_TYPE>& matrix,const Vec& rhs,Vec& result,bool syncPrecon) {
        setMatrix(matrix,syncPrecon);
        return solve(rhs,result);
    }
    SOLVER_RESULT solveMatrix(const FixedSparseMatrix<T,KERNEL_TYPE>& matrix,const Vec& rhs,Vec& result,bool syncPrecon) {
        setMatrix(matrix,syncPrecon);
        return solve(rhs,result);
    }
    virtual sizeType getIterationsCount() const {
        return _iterationsOut;
    }
    virtual T getResidual() const {
        return _residualOut;
    }
    virtual void setCallback(typename boost::shared_ptr<Callback<T,KERNEL_TYPE> > cb) {
        _cb=cb;
    }
protected:
    sizeType _iterationsOut;
    T _residualOut;
    boost::shared_ptr<Callback<T,KERNEL_TYPE> > _cb;
};

//============================================================================
// No Preconditioner Verbose
template <typename T,typename KERNEL_TYPE=Kernel<T> >
struct NoPreconSolver : public Solver<T,KERNEL_TYPE> {
public:
    typedef typename KERNEL_TYPE::Vec Vec;
    using Solver<T,KERNEL_TYPE>::_residualOut;
    using Solver<T,KERNEL_TYPE>::_iterationsOut;
    virtual void setMatrix(const FixedSparseMatrix<T,KERNEL_TYPE>& matrix,bool syncPrecon) {}
    virtual typename Solver<T, KERNEL_TYPE>::SOLVER_RESULT solve(const Vec& rhs,Vec& result) {
        KERNEL_TYPE::copy(rhs,result);

        _residualOut=0.0f;
        _iterationsOut=0;
        return Solver<T, KERNEL_TYPE>::SUCCESSFUL;
    }
};

//============================================================================
// Diagonal Preconditioner
template <typename T,typename KERNEL_TYPE=Kernel<T> >
struct DiagonalPreconSolver : public Solver<T,KERNEL_TYPE> {
public:
    typedef typename KERNEL_TYPE::Vec Vec;
    using Solver<T,KERNEL_TYPE>::_residualOut;
    using Solver<T,KERNEL_TYPE>::_iterationsOut;
    virtual void setMatrix(const FixedSparseMatrix<T,KERNEL_TYPE>& matrix,bool syncPrecon) {}
    virtual typename Solver<T, KERNEL_TYPE>::SOLVER_RESULT solve(const Vec& rhs,Vec& result) {
        KERNEL_TYPE::cwiseProd(rhs,_diag,result);

        _residualOut=0.0f;
        _iterationsOut=0;
        return Solver<T, KERNEL_TYPE>::SUCCESSFUL;
    }
	Vec _diag;
};

//============================================================================
// Diagonal Preconditioner
template <typename T,typename KERNEL_TYPE=Kernel<T> >
struct DirectPreconSolver : public Solver<T,KERNEL_TYPE> {
public:
    typedef typename KERNEL_TYPE::Vec Vec;
    using Solver<T,KERNEL_TYPE>::_residualOut;
    using Solver<T,KERNEL_TYPE>::_iterationsOut;
	virtual void setMatrix(const FixedSparseMatrix<T,KERNEL_TYPE>& matrix,bool syncPrecon) {
		Eigen::SparseMatrix<T,0,sizeType> emat;
		matrix.toEigen(emat);
		setMatrix(emat,syncPrecon);
	}
	virtual void setMatrix(const Eigen::SparseMatrix<T,0,sizeType>& matrix,bool syncPrecon) {_fact.compute(matrix);}
    virtual typename Solver<T, KERNEL_TYPE>::SOLVER_RESULT solve(const Vec& rhs,Vec& result) {
        result=_fact.solve(rhs);

        _residualOut=0.0f;
        _iterationsOut=0;
        return Solver<T, KERNEL_TYPE>::SUCCESSFUL;
    }
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<T,0,sizeType> > _fact;
};

//============================================================================
// If M is known, just solve by multiplication: Mb=x
template <typename T,typename KERNEL_TYPE=Kernel<T> >
struct ExplicitPreconSolver : public Solver<T,KERNEL_TYPE> {
public:
    typedef typename KERNEL_TYPE::Vec Vec;
    using Solver<T,KERNEL_TYPE>::_residualOut;
    using Solver<T,KERNEL_TYPE>::_iterationsOut;
    virtual void setMatrix(const FixedSparseMatrix<T,KERNEL_TYPE>& matrix,bool syncPrecon) {
        _matrix.reset(new DefaultKrylovMatrix<T,KERNEL_TYPE>(matrix));
    }
    virtual void setKrylovMatrix(boost::shared_ptr<KrylovMatrix<T,KERNEL_TYPE> > matrix) {
        _matrix=matrix;
    }
    typename Solver<T, KERNEL_TYPE>::SOLVER_RESULT solve(const Vec& rhs,Vec& result) {
        _matrix->multiply(rhs,result);
        return Solver<T, KERNEL_TYPE>::SUCCESSFUL;
    }
protected:
    boost::shared_ptr<KrylovMatrix<T,KERNEL_TYPE> > _matrix;
};

//============================================================================
// Incomplete Cholesky factorization, level zero, with option for modified version.
// Set modification_parameter between zero (regular incomplete Cholesky) and
// one (fully modified version), with values close to one usually giving the best
// results. The min_diagonal_ratio parameter is used to detect and correct
// problems in factorization: if a pivot is this much less than the diagonal
// entry from the original matrix, the original matrix entry is used instead.
template <typename T,typename KERNEL_TYPE=Kernel<T> >
struct MIC0Solver : public Solver<T,KERNEL_TYPE> {
public:
    template <typename T2>
    struct SparseColumnLowerFactor {
    public:
        typedef typename KERNEL_TYPE::template Rebind<T2> R2;
        typedef typename R2::value::Vec VecT2;
        //using Solver<T,KERNEL_TYPE>::_residualOut;
        //using Solver<T,KERNEL_TYPE>::_iterationsOut;
        SparseColumnLowerFactor(sizeType n=0)
            :_n(n),_invDiag(n),_colStart(n+1),_adiag(n) {}
        void clear(void) {
            _n=0;
            _invDiag.clear();
            _value.clear();
            _rowIndex.clear();
            _colStart.clear();
            _adiag.clear();
        }
        void resize(sizeType n) {
            _n=n;
            _invDiag.resize(n);
            _colStart.resize(n+1);
            _adiag.resize(n);
        }
        void writeMatlab(std::ostream &output,const string& vName) {
            output << vName << "=sparse([";
            for(sizeType i=0; i<_n; ++i) {
                output << " " << i+1;
                for(sizeType j=_colStart[i]; j<_colStart[i+1]; ++j) {
                    output << " " << _rowIndex[j]+1;
                }
            }

            output << "],...\n  [";
            for(sizeType i=0; i<_n; ++i) {
                output << " " << i+1;
                for(sizeType j=_colStart[i]; j<_colStart[i+1]; ++j) {
                    output << " " << i+1;
                }
            }

            output << "],...\n  [";
            for(sizeType i=0; i<_n; ++i) {
                output << " " << (_invDiag[i]!=0.0f ? 1.0f/_invDiag[i] : 0.0f);
                for(sizeType j=_colStart[i]; j<_colStart[i+1]; ++j) {
                    output << " " << _value[j];
                }
            }

            output << "], " << _n << ", " << _n << ");" <<std::endl;
        }
    public:
        sizeType _n;
        VecT2 _invDiag; // reciprocals of diagonal elements
        std::vector<T2> _value; // values below the diagonal, listed column by column
        std::vector<sizeType> _rowIndex; // a list of all row indices, for each column in turn
        std::vector<sizeType> _colStart; // where each column begins in rowindex (plus an extra entry at the end, of #nonzeros)
        VecT2 _adiag; // just used in factorization: minimum "safe" diagonal entry allowed
    };
    struct MatrixFilter
    {
        virtual bool operator()(sizeType r,sizeType c) const{return false;}
    };
public:
    typedef typename KERNEL_TYPE::Vec Vec;
    using Solver<T,KERNEL_TYPE>::_residualOut;
    using Solver<T,KERNEL_TYPE>::_iterationsOut;
    virtual void setMatrix(const FixedSparseMatrix<T,KERNEL_TYPE>& matrix,bool syncPrecon) {
        MIC0(matrix,_icFactor);
    }
    virtual void setFilter(boost::shared_ptr<MatrixFilter> filter){_filter=filter;}
    virtual typename Solver<T,KERNEL_TYPE>::SOLVER_RESULT solve(const Vec& rhs,Vec& result) {
        solveLower(_icFactor,rhs,result);
        solveLowerTransposeInPlace(_icFactor,result);

        _residualOut=0.0f;
        _iterationsOut=1;
        return Solver<T,KERNEL_TYPE>::SUCCESSFUL;
    }
protected:
    void MIC0(const FixedSparseMatrix<T,KERNEL_TYPE> &matrix,
              SparseColumnLowerFactor<T> &factor,
              T modificationParameter=0.97f,
              T minDiagonalRatio=0.25f) {
        // first copy lower triangle of matrix into factor (Note: assuming A is symmetric of course!)
        factor.resize(matrix.rows());
        KERNEL_TYPE::zero(factor._invDiag); // important: eliminate old values from previous solves!
        factor._value.resize(0);
        factor._rowIndex.resize(0);
        KERNEL_TYPE::zero(factor._adiag);

        for(sizeType i=0; i<matrix.rows(); ++i) {
            factor._colStart[i]=(sizeType)factor._rowIndex.size();
            ConstSMIterator<T> beg=matrix.begin(i),end=matrix.end(i);
            while(beg!=end && _filter && (*_filter)(beg.row(),beg.col()))
                ++beg;
            for(;beg!=end;) {
                if(beg.col()>i) {
                    factor._rowIndex.push_back(beg.col());
                    factor._value.push_back(*beg);
                } else if(beg.col()==i) {
                    factor._invDiag[i]=factor._adiag[i]=*beg;
                }
                //++beg;
                do{
                    ++beg;
                }while(beg!=end && _filter && (*_filter)(beg.row(),beg.col()));
            }
        }
        factor._colStart[matrix.rows()]=(sizeType)factor._rowIndex.size();
        // now do the incomplete factorization (figure out numerical values)

        // MATLAB code:
        // L=tril(A);
        // for k=1:size(L,2)
        //   L(k,k)=sqrt(L(k,k));
        //   L(k+1:end,k)=L(k+1:end,k)/L(k,k);
        //   for j=find(L(:,k))'
        //     if j>k
        //       fullupdate=L(:,k)*L(j,k);
        //       incompleteupdate=fullupdate.*(A(:,j)~=0);
        //       missing=sum(fullupdate-incompleteupdate);
        //       L(j:end,j)=L(j:end,j)-incompleteupdate(j:end);
        //       L(j,j)=L(j,j)-omega*missing;
        //     end
        //   end
        // end

        for(sizeType k=0; k<matrix.rows(); ++k) {
            if(factor._adiag[k]==0) continue; // null row/column

            // figure out the final L(k,k) entry
            if(factor._invDiag[k]<minDiagonalRatio*factor._adiag[k])
                factor._invDiag[k]=1/sqrt(factor._adiag[k]); // drop to Gauss-Seidel here if the pivot looks dangerously small
            else
                factor._invDiag[k]=1/sqrt(factor._invDiag[k]);

            // finalize the k'th column L(:,k)
            for(sizeType p=factor._colStart[k]; p<factor._colStart[k+1]; ++p) {
                factor._value[p]*=factor._invDiag[k];
            }

            // incompletely eliminate L(:,k) from future columns, modifying diagonals
            for(sizeType p=factor._colStart[k]; p<factor._colStart[k+1]; ++p) {
                sizeType j=factor._rowIndex[p]; // work on column j
                T multiplier=factor._value[p];
                T missing=0;
                sizeType a=factor._colStart[k];

                // first look for contributions to missing from dropped entries above the diagonal in column j
                sizeType b=0;
                ConstSMIterator<T> beg=matrix.begin(j),end=matrix.end(j);
                while(beg!=end && _filter && (*_filter)(beg.row(),beg.col()))
                    ++beg;
                while(a<factor._colStart[k+1] && factor._rowIndex[a]<j) {
                    // look for factor.rowindex[a] in matrix.index[j] starting at b
                    while(beg!=end) {
                        if(beg.col()<factor._rowIndex[a]) {
                            do{
                                ++beg;
                            }while(beg!=end && _filter && (*_filter)(beg.row(),beg.col()));
                            //++beg;
                            ++b;
                        } else if(beg.col()==factor._rowIndex[a])
                            break;
                        else {
                            missing+=factor._value[a];
                            break;
                        }
                    }
                    ++a;
                }

                // adjust the diagonal j,j entry
                if(a<factor._colStart[k+1] && factor._rowIndex[a]==j) {
                    factor._invDiag[j]-=multiplier*factor._value[a];
                }
                ++a;

                // and now eliminate from the nonzero entries below the diagonal in column j (or add to missing if we can't)
                b=factor._colStart[j];
                while(a<factor._colStart[k+1] && b<factor._colStart[j+1]) {
                    if(factor._rowIndex[b]<factor._rowIndex[a])
                        ++b;
                    else if(factor._rowIndex[b]==factor._rowIndex[a]) {
                        factor._value[b]-=multiplier*factor._value[a];
                        ++a;
                        ++b;
                    } else {
                        missing+=factor._value[a];
                        ++a;
                    }
                }

                // and if there's anything left to do, add it to missing
                while(a<factor._colStart[k+1]) {
                    missing+=factor._value[a];
                    ++a;
                }

                // and do the final diagonal adjustment from the missing entries
                factor._invDiag[j]-=modificationParameter*multiplier*missing;
            }
        }
    }
    // Solution routines with lower triangular matrix.
    // solve L*result=rhs
    void solveLower(const SparseColumnLowerFactor<T> &factor,const Vec& rhs,Vec& result) {
        ASSERT(factor._n == rhs.size());
        ASSERT(factor._n == result.size());

        result=rhs;
        for(sizeType i=0; i<factor._n; ++i) {
            result[i]*=factor._invDiag[i];
            for(sizeType j=factor._colStart[i]; j<factor._colStart[i+1]; ++j) {
                result[factor._rowIndex[j]]-=factor._value[j]*result[i];
            }
        }
    }
    // solve L^T*result=rhs
    void solveLowerTransposeInPlace(const SparseColumnLowerFactor<T> &factor,Vec& x) {
        ASSERT(factor._n == x.size());
        ASSERT(factor._n > 0);

        sizeType i=factor._n;
        do {
            --i;
            for(sizeType j=factor._colStart[i]; j<factor._colStart[i+1]; ++j) {
                x[i]-=factor._value[j]*x[factor._rowIndex[j]];
            }
            x[i]*=factor._invDiag[i];
        } while(i!=0);
    }
    SparseColumnLowerFactor<T> _icFactor;
    boost::shared_ptr<MatrixFilter> _filter;
};

//============================================================================
// Encapsulates the Conjugate Gradient algorithm with incomplete Cholesky
// factorization preconditioner.
template <typename KERNEL_TYPE>
struct Nullspace {
public:
    typedef typename KERNEL_TYPE::Vec Vec;
    virtual void projectOut(Vec& r) =0;
};
template <typename T,typename KERNEL_TYPE=Kernel<T>,typename PRECONDITIONER_TYPE=MIC0Solver<T> >
struct PCGSolver : public Solver<T,KERNEL_TYPE> {
public:
    typedef typename KERNEL_TYPE::Vec Vec;
    using Solver<T,KERNEL_TYPE>::_cb;
    using Solver<T,KERNEL_TYPE>::_residualOut;
    using Solver<T,KERNEL_TYPE>::_iterationsOut;
    PCGSolver(void) {
        setSolverParameters(1e-5f,1000);
        _ns=NULL;
        _useIPCG=false;
        _oneNorm=false;
        _warmStart=false;
    }
    void setNullspace(Nullspace<KERNEL_TYPE>* ns) {
        _ns=ns;
    }
    virtual void setSolverParameters(T toleranceFactor,sizeType maxIterations) {
        _toleranceFactor=toleranceFactor;
        if(_toleranceFactor<1e-30f)
            _toleranceFactor=1e-30f;
        _maxIterations=maxIterations;
    }
    virtual void resize(const sizeType& n) {
        if(_r.size()!=n) {
            _s.resize(n);
            _z.resize(n);
            _r.resize(n);
            if(_useIPCG)
                _m.resize(n);
        }
    }
    virtual void setMatrix(const FixedSparseMatrix<T,KERNEL_TYPE>& matrix,bool syncPrecon) {
        resize(matrix.rows());
        if(syncPrecon)_pre.setMatrix(matrix,true);
        _fixedMatrix.reset(new DefaultKrylovMatrix<T,KERNEL_TYPE>(matrix));
    }
    virtual void setKrylovMatrix(boost::shared_ptr<KrylovMatrix<T,KERNEL_TYPE> > matrix) {
        resize(matrix->n());
        _fixedMatrix=matrix;
    }
    virtual typename Solver<T,KERNEL_TYPE>::SOLVER_RESULT solve(const Vec& rhs,Vec& result) {
        if(_cb)
            _cb->reset();

        if(_warmStart){
            _fixedMatrix->multiply(result,_r);
            KERNEL_TYPE::sub(rhs,_r,_r);
        }else{
            result.resize(_fixedMatrix->n());
            KERNEL_TYPE::zero(result);
            KERNEL_TYPE::copy(rhs,_r);
        }
        if(_ns)
            _ns->projectOut(_r);
        _residualOut=_oneNorm ? KERNEL_TYPE::absMax(_r) : KERNEL_TYPE::norm(_r);
        if(_residualOut==0.0f) {
            _iterationsOut=0;
            return Solver<T,KERNEL_TYPE>::SUCCESSFUL;
        }
        _pre.solve(_r,_z);
        T tol=_toleranceFactor*_residualOut;

        T rho=KERNEL_TYPE::dot(_z,_r);
        if(rho==0 || rho!=rho) {
            _iterationsOut=0;
            return Solver<T,KERNEL_TYPE>::NOT_CONVERGENT;
        }

        KERNEL_TYPE::copy(_z,_s);
        sizeType iteration;
        for(iteration=0; iteration<_maxIterations; ++iteration) {
            _fixedMatrix->multiply(_s,_z);
            T alpha=rho/KERNEL_TYPE::dot(_s,_z);
            KERNEL_TYPE::addScaled(alpha,_s,result);
            if(_useIPCG)
                KERNEL_TYPE::copy(_r,_m);
            KERNEL_TYPE::addScaled(-alpha,_z,_r);

            if(_ns)
                _ns->projectOut(_r);

            _residualOut=_oneNorm ? KERNEL_TYPE::absMax(_r) : KERNEL_TYPE::norm(_r);
            if(_cb){
                if((*_cb)(result,_residualOut,iteration))
                    return Solver<T,KERNEL_TYPE>::USER_REQUEST;
            }
            if(_residualOut<=tol) {
                _iterationsOut=iteration+1;
                return Solver<T,KERNEL_TYPE>::SUCCESSFUL;
            }
            if(_pre.solve(_r,_z) != Solver<T,KERNEL_TYPE>::SUCCESSFUL)
                return Solver<T,KERNEL_TYPE>::NOT_CONVERGENT;

            T rhoNew=KERNEL_TYPE::dot(_z,_r);
            T beta;
            if(_useIPCG)
                beta=(rhoNew-KERNEL_TYPE::dot(_z,_m))/rho;
            else beta=rhoNew/rho;
            KERNEL_TYPE::addScaled(beta,_s,_z);
            KERNEL_TYPE::copy(_z,_s); // s=beta*s+z
            rho=rhoNew;
        }

        _iterationsOut=iteration;
        return Solver<T,KERNEL_TYPE>::NOT_CONVERGENT;
    }
    virtual void setUseIPCG(bool useIPCG){_useIPCG=useIPCG;}
    Solver<T,KERNEL_TYPE>* getPre() {return &_pre;}
    const Solver<T,KERNEL_TYPE>* getPre() const {return &_pre;}
    void setUseOneNorm(bool oneNorm){_oneNorm=oneNorm;}
    void setUseWarmStart(bool warmStart){_warmStart=warmStart;}
protected:
    // internal structures
    PRECONDITIONER_TYPE _pre;
    //_m is for flexible CG
    Vec _z, _s, _r, _m; // temporary vectors for PCG
    boost::shared_ptr<KrylovMatrix<T,KERNEL_TYPE> > _fixedMatrix;
    // parameters
    T _toleranceFactor;
    sizeType _maxIterations;
    //null space
    Nullspace<KERNEL_TYPE>* _ns;
    //do you want to use complicated preconditioner such MG?
    bool _useIPCG;
    bool _oneNorm;
    bool _warmStart;
};

//============================================================================
// Encapsulates the BiConjugate Gradient Stabilized algorithm with incomplete Cholesky
// factorization preconditioner.
/*
//============================================================================
//reference code
template < class Matrix, class Vector, class Preconditioner, class Real >
int
BiCGSTAB(const Matrix &A, Vector &x, const Vector &b,
         const Preconditioner &M, int &max_iter, Real &tol)
{
  Real resid;
  Vector rho_1(1), rho_2(1), alpha(1), beta(1), omega(1);
  Vector p, phat, s, shat, t, v;

  Real normb = norm(b);
  Vector r = b - A * x;
  Vector rtilde = r;

  if (normb == 0.0)
    normb = 1;

  if ((resid = norm(r) / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  for (int i = 1; i <= max_iter; i++) {
    rho_1(0) = dot(rtilde, r);
    if (rho_1(0) == 0) {
      tol = norm(r) / normb;
      return 2;
    }
    if (i == 1)
      p = r;
    else {
      beta(0) = (rho_1(0)/rho_2(0)) * (alpha(0)/omega(0));
      p = r + beta(0) * (p - omega(0) * v);
    }
    phat = M.solve(p);
    v = A * phat;
    alpha(0) = rho_1(0) / dot(rtilde, v);
    s = r - alpha(0) * v;
    if ((resid = norm(s)/normb) < tol) {
      x += alpha(0) * phat;
      tol = resid;
      return 0;
    }
    shat = M.solve(s);
    t = A * shat;
    omega = dot(t,s) / dot(t,t);
    x += alpha(0) * phat + omega(0) * shat;
    r = s - omega(0) * t;

    rho_2(0) = rho_1(0);
    if ((resid = norm(r) / normb) < tol) {
      tol = resid;
      max_iter = i;
      return 0;
    }
    if (omega(0) == 0) {
      tol = norm(r) / normb;
      return 3;
    }
  }

  tol = resid;
  return 1;
}
*/
template <typename T,typename KERNEL_TYPE=Kernel<T>,typename PRECONDITIONER_TYPE=MIC0Solver<T> >
struct BiCGStabSolver : public Solver<T,KERNEL_TYPE> {
public:
    typedef typename KERNEL_TYPE::Vec Vec;
    using Solver<T,KERNEL_TYPE>::_cb;
    using Solver<T,KERNEL_TYPE>::_residualOut;
    using Solver<T,KERNEL_TYPE>::_iterationsOut;
    BiCGStabSolver(void) {
        setSolverParameters(1e-5f,1000);
        _oneNorm=true;
    }
    virtual void setSolverParameters(T toleranceFactor,sizeType maxIterations) {
        _toleranceFactor=toleranceFactor;
        if(_toleranceFactor<1e-30f)
            _toleranceFactor=1e-30f;
        _maxIterations=maxIterations;
    }
    virtual void resize(const sizeType& n) {
        if(_r.size()!=n) {
            _r.resize(n);
            _rHat.resize(n);
            _v.resize(n);
            _p.resize(n);
            _pHat.resize(n);
            _s.resize(n);
            _sHat.resize(n);
            _t.resize(n);
        }
    }
    virtual void setMatrix(const FixedSparseMatrix<T,KERNEL_TYPE> &matrix,bool syncPrecon) {
        resize(matrix.rows());
        if(syncPrecon)_pre.setMatrix(matrix,true);
        _fixedMatrix.reset(new DefaultKrylovMatrix<T,KERNEL_TYPE>(matrix));
    }
    virtual void setKrylovMatrix(boost::shared_ptr<KrylovMatrix<T,KERNEL_TYPE> > matrix) {
        resize(matrix->n());
        _fixedMatrix=matrix;
    }
    virtual typename Solver<T,KERNEL_TYPE>::SOLVER_RESULT solve(const Vec& rhs,Vec& result) {
        if(_cb)
            _cb->reset();

        result.resize(_fixedMatrix->n());
        _fixedMatrix->multiply(result,_r);
        KERNEL_TYPE::sub(rhs,_r,_r);
        KERNEL_TYPE::copy(_r,_rHat);
        _residualOut=_oneNorm ? KERNEL_TYPE::absMax(_r) : KERNEL_TYPE::norm(_r);
        if(_residualOut==0.0f) {
            _iterationsOut=0;
            return Solver<T,KERNEL_TYPE>::SUCCESSFUL;
        }
        T tol=_toleranceFactor*_residualOut;

        T rho=1.0f;
        T alpha=1.0f;
        T omega=1.0f;
        T beta;

        sizeType iteration;
        for(iteration=0; iteration<_maxIterations; ++iteration) {
            T newRho=KERNEL_TYPE::dot(_rHat,_r);
            beta=(newRho/rho)*(alpha/omega);
            rho=newRho;	//rho_i-1 -> rho_i

            //_p=_r+(_p-_v*omega)*beta;	//p_i-1 -> p_i
            if(iteration == 0) {
                _p=_r;
            } else {
                KERNEL_TYPE::addScaled(-omega,_v,_p);
                KERNEL_TYPE::scale(beta,_p);
                KERNEL_TYPE::addScaled(1.0f,_r,_p);
            }

            if(_pre.solve(_p,_pHat) != Solver<T,KERNEL_TYPE>::SUCCESSFUL)
                return Solver<T,KERNEL_TYPE>::NOT_CONVERGENT;
            _fixedMatrix->multiply(_pHat,_v);	//v_i-1 -> v_i
            alpha=rho/KERNEL_TYPE::dot(_rHat,_v);

            //_s=_r-_v*alpha;
            KERNEL_TYPE::copy(_r,_s);
            KERNEL_TYPE::addScaled(-alpha,_v,_s);

            if(_pre.solve(_s,_sHat) != Solver<T,KERNEL_TYPE>::SUCCESSFUL)
                return Solver<T,KERNEL_TYPE>::NOT_CONVERGENT;
            _fixedMatrix->multiply(_sHat,_t);
            omega=KERNEL_TYPE::dot(_t,_s)/
                  KERNEL_TYPE::dot(_t,_t);	//omega_i-1 -> omega

            //_x=_x+_y*alpha+_z*omega;	//x_i-1 -> x_i
            KERNEL_TYPE::addScaled(alpha,_pHat,result);
            KERNEL_TYPE::addScaled(omega,_sHat,result);

            //_r=_s-_t*omega;	//r_i-1 -> r_i
            KERNEL_TYPE::copy(_s,_r);
            KERNEL_TYPE::addScaled(-omega,_t,_r);
            
            _residualOut=_oneNorm ? KERNEL_TYPE::absMax(_r) : KERNEL_TYPE::norm(_r);
            if(_cb){
                if((*_cb)(result,_residualOut,iteration))
                    return Solver<T,KERNEL_TYPE>::USER_REQUEST;
            }
            if(_residualOut<=tol) {
                _iterationsOut=iteration+1;
                return Solver<T,KERNEL_TYPE>::SUCCESSFUL;
            }
        }

        _iterationsOut=iteration;
        return Solver<T,KERNEL_TYPE>::NOT_CONVERGENT;
    }
    Solver<T,KERNEL_TYPE>* getPre() {return &_pre;}
    const Solver<T,KERNEL_TYPE>* getPre() const {return &_pre;}
    void setUseOneNorm(bool oneNorm){_oneNorm=oneNorm;}
protected:
    // internal structures
    PRECONDITIONER_TYPE _pre;
    Vec _r, _rHat, _v, _p, _pHat, _s, _sHat, _t; // temporary vectors for PCG
    boost::shared_ptr<KrylovMatrix<T,KERNEL_TYPE> > _fixedMatrix;
    // parameters
    T _toleranceFactor;
    sizeType _maxIterations;
    bool _oneNorm;
};

//============================================================================
// Encapsulates the GMRES algorithm, just a port of the code tminres from IML
template <typename T,typename KERNEL_TYPE=Kernel<T>,typename PRECONDITIONER_TYPE=NoPreconSolver<T> >
struct GMRESSolver : public Solver<T,KERNEL_TYPE> {
public:
    typedef typename KERNEL_TYPE::Vec Vec;
    using Solver<T,KERNEL_TYPE>::_cb;
    using Solver<T,KERNEL_TYPE>::_residualOut;
    using Solver<T,KERNEL_TYPE>::_iterationsOut;
    GMRESSolver(void) {
        setSolverParameters(1e-5f,1000);
    }
    virtual void setSolverParameters(T toleranceFactor,sizeType maxIterations) {
        _toleranceFactor=toleranceFactor;
        if(_toleranceFactor<1e-30f)
            _toleranceFactor=1e-30f;
        _maxIterations=maxIterations;
        _m=40;
    }
    virtual void resize(const sizeType& n) {
        if(_r.size()!=n) {
            _s.resize(_m+1);
            _cs.resize(_m+1);
            _sn.resize(_m+1);
            _w.resize(n);
            _r.resize(n);
            _tmp.resize(n);
            _tmp2.resize(n);
            _v.resize(_m+1);
            for(sizeType i=0; i<(sizeType)_v.size(); i++)
                _v[i].resize(n);
            _H.resize(_m+1,_m+1);
        }
    }
    virtual void setMatrix(const FixedSparseMatrix<T,KERNEL_TYPE> &matrix,bool syncPrecon) {
        resize(matrix.rows());
        if(syncPrecon)_pre.setMatrix(matrix,true);
        _fixedMatrix.reset(new DefaultKrylovMatrix<T,KERNEL_TYPE>(matrix));
    }
    virtual void setKrylovMatrix(typename boost::shared_ptr<KrylovMatrix<T,KERNEL_TYPE> > matrix) {
        resize(matrix->n());
        _fixedMatrix=matrix;
    }
    virtual typename Solver<T,KERNEL_TYPE>::SOLVER_RESULT solve(const Vec& rhs,Vec& result) {
        if(_cb)
            _cb->reset();

        sizeType i,j=1,k;
        result.resize(_fixedMatrix->n());
        KERNEL_TYPE::zero(result);
        _pre.solve(rhs,_tmp);
        T normb=_tmp.norm();

        _fixedMatrix->multiply(result,_tmp);
        KERNEL_TYPE::sub(rhs,_tmp,_tmp2);
        _pre.solve(_tmp2,_r);
        T beta=_r.norm();

        if(normb == 0.0f)
            normb=1;

        if((_residualOut=beta/normb) <= _toleranceFactor) {
            _iterationsOut=0;
            return Solver<T,KERNEL_TYPE>::SUCCESSFUL;
        }

        while(j <= _maxIterations) {
            _v[0]=_r*(1.0f/beta);    // ??? r / beta
            KERNEL_TYPE::zero(_s);
            _s(0)=beta;

            if(_cb){
                if((*_cb)(result,_residualOut,j))
                    return Solver<T,KERNEL_TYPE>::USER_REQUEST;
            }
            for(i=0; i < _m && j <= _maxIterations; i++,j++) {
                _fixedMatrix->multiply(_v[i],_tmp);
                _pre.solve(_tmp,_w);
                for(k=0; k <= i; k++) {
                    _H(k,i)=KERNEL_TYPE::dot(_w,_v[k]);
                    _w-=_H(k,i)*_v[k];
                }
                _H(i+1,i)=_w.norm();
                _v[i+1]=_w*(1.0f/_H(i+1,i)); // ??? w / H(i+1, i)

                for(k=0; k < i; k++)
                    applyPlaneRotation(_H(k,i),_H(k+1,i),_cs(k),_sn(k));

                generatePlaneRotation(_H(i,i),_H(i+1,i),_cs(i),_sn(i));
                applyPlaneRotation(_H(i,i),_H(i+1,i),_cs(i),_sn(i));
                applyPlaneRotation(_s(i),_s(i+1),_cs(i),_sn(i));

                if((_residualOut=std::abs(_s(i+1))/normb) < _toleranceFactor) {
                    update(result,i,_H,_s,_v);
                    _iterationsOut=j;
                    return Solver<T,KERNEL_TYPE>::SUCCESSFUL;
                }
            }

            update(result,_m-1,_H,_s,_v);
            _fixedMatrix->multiply(result,_tmp);
            KERNEL_TYPE::sub(rhs,_tmp,_tmp2);
            _pre.solve(_tmp2,_r);
            beta=_r.norm();
            if((_residualOut=beta/normb) < _toleranceFactor) {
                _iterationsOut=j;
                return Solver<T,KERNEL_TYPE>::SUCCESSFUL;
            }
        }

        return Solver<T,KERNEL_TYPE>::NOT_CONVERGENT;
    }
    Solver<T,KERNEL_TYPE>* getPre() {return &_pre;}
    const Solver<T,KERNEL_TYPE>* getPre() const {return &_pre;}
protected:
    void update(Vec& x,const sizeType& k,const Eigen::Matrix<T,-1,-1> &h,const Vec &s,
                const std::vector<Vec,Eigen::aligned_allocator<Vec> >& v) const {
        Vec y;
        KERNEL_TYPE::copy(s,y);
        // Backsolve:
        for(sizeType i=k; i>=0; i--) {
            y(i)/=h(i,i);
            for(sizeType j=i-1; j>=0; j--)
                y(j)-=h(j,i)*y(i);
        }
        for(sizeType j=0; j<=k; j++)
            x+=v[j]*y(j);
    }
    void generatePlaneRotation(const T &dx,const T &dy,T &cs,T &sn) const {
        if(dy == 0.0f) {
            cs=1.0f;
            sn=0.0f;
        } else if(std::abs(dy) > std::abs(dx)) {
            T temp=dx/dy;
            sn=1.0f/sqrt(1.0f+temp*temp);
            cs=temp*sn;
        } else {
            T temp=dy/dx;
            cs=1.0f/sqrt(1.0f+temp*temp);
            sn=temp*cs;
        }
    }
    void applyPlaneRotation(T &dx,T &dy,const T &cs,const T &sn) const {
        T temp=cs*dx+sn*dy;
        dy=-sn*dx+cs*dy;
        dx=temp;
    }
    //internal structure
    PRECONDITIONER_TYPE _pre;
    Vec _s,_cs,_sn,_w,_r,_tmp,_tmp2;
    std::vector<Vec,Eigen::aligned_allocator<Vec> > _v;
    Eigen::Matrix<T,-1,-1> _H;
    boost::shared_ptr<KrylovMatrix<T,KERNEL_TYPE> > _fixedMatrix;
    // parameters
    T _toleranceFactor;
    sizeType _maxIterations;
    sizeType _m;
};

PRJ_END

#endif
