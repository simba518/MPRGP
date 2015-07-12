#ifndef DACG_SOLVER_H
#define DACG_SOLVER_H

#include "LinearSolver.h"

PRJ_BEGIN

template <typename T,typename KERNEL_TYPE=Kernel<T> >
struct ExactInvSolver : public Solver<T,KERNEL_TYPE> {
public:
    typedef typename KERNEL_TYPE::Vec Vec;
    using Solver<T,KERNEL_TYPE>::_residualOut;
    using Solver<T,KERNEL_TYPE>::_iterationsOut;
    virtual void setMatrix(const FixedSparseMatrix<T,KERNEL_TYPE>& matrix,bool syncPrecon) {
		Eigen::SparseMatrix<T,0,sizeType> m;
		matrix.toEigen(m);
		_sol.compute(m);
		_dense=false;
	}
    virtual void setMatrix(const Eigen::Matrix<T,-1,-1>& matrix,bool syncPrecon) {
		_invM=matrix.inverse();
		_dense=true;
	}
    virtual typename Solver<T, KERNEL_TYPE>::SOLVER_RESULT solve(const Vec& rhs,Vec& result) {
		if(_dense)result=_invM*rhs;
		else result=_sol.solve(rhs);
        _residualOut=0.0f;
        _iterationsOut=0;
        return Solver<T, KERNEL_TYPE>::SUCCESSFUL;
    }
protected:
	Eigen::SimplicialCholesky<Eigen::SparseMatrix<T,0,sizeType> > _sol;
	Eigen::Matrix<T,-1,-1> _invM;
	bool _dense;
};
template <typename T,typename KERNEL_TYPE=Kernel<T>,typename PRECONDITIONER_TYPE=NoPreconSolver<T> >
class DACGSolver
{
public:
	typedef typename KERNEL_TYPE::Vec Vec;
	typedef typename Solver<T,KERNEL_TYPE>::SOLVER_RESULT RESULT;
	//setter
	DACGSolver(void){setSolverParameters(1e-20f,1e-5f,10000);}
	virtual void setSolverParameters(T eps1,T eps2,sizeType maxIterations)
	{
		_eps1=std::min<T>(std::max<T>(eps1,0.0f),1.0f);
		_eps2=std::min<T>(std::max<T>(eps2,0.0f),1.0f);
		_maxIter=std::max<sizeType>(maxIterations,100);
		_maxIterInner=1000;
	}
	virtual Solver<T,KERNEL_TYPE>* getPre() {return &_pre;}
    virtual const Solver<T,KERNEL_TYPE>* getPre() const {return &_pre;}
	virtual void setCallback(typename boost::shared_ptr<Callback<T,KERNEL_TYPE> > cb) {_cb=cb;}
	virtual void resize(const sizeType& n) 
	{
		if(_xk.size() != n)
		{
			_xk.resize(n);
			_xkA.resize(n);
			_xkB.resize(n);

			_pk.resize(n);
			_pkA.resize(n);
			_pkB.resize(n);
			_pkTilde.resize(n);

			_gk.resize(n);
			_gkM.resize(n);
		}
    }
	//set matrix A
	virtual void setA(const FixedSparseMatrix<T,KERNEL_TYPE>& A,bool syncPrecon)
	{
		resize(A.rows());
		_A.reset(new DefaultKrylovMatrix<T,KERNEL_TYPE>(A));
		if(syncPrecon)_pre.setMatrix(A,syncPrecon);
	}
    virtual void setKrylovA(boost::shared_ptr<KrylovMatrix<T,KERNEL_TYPE> > A){resize(A->n());_A=A;}
	virtual void setKrylovA(const Eigen::SparseMatrix<T,0,sizeType>& A)
	{
		resize(A.rows());
		_A.reset(new DefaultEigenKrylovMatrix<T,KERNEL_TYPE>(A));
	}
	virtual void setKrylovA(const Eigen::Matrix<T,-1,-1>& A)
	{
		resize(A.rows());
		_A.reset(new DefaultDenseEigenKrylovMatrix<T,KERNEL_TYPE>(A));
	}
	//set matrix B
	virtual void setB(const FixedSparseMatrix<T,KERNEL_TYPE>& B)
	{
		resize(B.rows());
		_B.reset(new DefaultKrylovMatrix<T,KERNEL_TYPE>(B));
	}
    virtual void setKrylovB(boost::shared_ptr<KrylovMatrix<T,KERNEL_TYPE> > B){resize(B->n());_B=B;}
	virtual void setKrylovB(const Eigen::SparseMatrix<T,0,sizeType>& B)
	{
		resize(B.rows());
		_B.reset(new DefaultEigenKrylovMatrix<T,KERNEL_TYPE>(B));
	}
	virtual void setKrylovB(const Eigen::Matrix<T,-1,-1>& B)
	{
		resize(B.rows());
		_B.reset(new DefaultDenseEigenKrylovMatrix<T,KERNEL_TYPE>(B));
	}
	//set U0
	virtual void setU0(const Eigen::Matrix<T,-1,-1>& U0)
	{
		_U0=U0;
		if(!_B)
		{
			setKrylovB(boost::shared_ptr<KrylovMatrix<T,KERNEL_TYPE> >(new IdentityKrylovMatrix<T,KERNEL_TYPE>(_A->n())));
		}
		//if(!_B)
		//{
		//	NOTIFY_MSG("U0 must be given after B!")
		//	exit(EXIT_FAILURE);
		//}
		Vec UC(_B->n()),UCB(_B->n());
		for(sizeType eig=1;eig<=_U0.cols();eig++)
		{
			INFOV("%dth Basis",eig)
			projectOut(_U0,eig,_U0.col(eig-1),UC,false);
			_B->multiply(UC,UCB);
			_U0.col(eig-1)=UC/std::sqrt(KERNEL_TYPE::dot(UC,UCB));
			//check orthogonality
			for(sizeType j=1;j<eig;j++)
			{
				_B->multiply(_U0.col(eig-1),UCB);
				INFOV("\tOrthogonal residual: %f",UCB.dot(_U0.col(j-1)))
			}
		}
	}
public:
	//solver
	void checkSolution(Vec& lambda,Eigen::Matrix<T,-1,-1>& U) const
	{
		checkMatrix();
		if(lambda.rows() != U.cols() || U.rows() != _A->n())
		{
			NOTIFY_MSG("Incorrect input size!")
			exit(EXIT_FAILURE);
		}
		Vec Ax,lambdaBx,err;
		Ax.resize(U.rows());
		lambdaBx.resize(U.rows());
		err.resize(U.rows());
		for(sizeType r=0;r<lambda.rows();r++)
		{
			_A->multiply(U.col(r),Ax);
			_B->multiply(U.col(r),lambdaBx);
			KERNEL_TYPE::scale(lambda[r],lambdaBx);
			KERNEL_TYPE::sub(Ax,lambdaBx,err);
			INFOV("lambda: %f, Ax: %f, lambdaBx: %f, err: %f",
				  lambda[r],
				  KERNEL_TYPE::norm(Ax),
				  KERNEL_TYPE::norm(lambdaBx),
				  KERNEL_TYPE::norm(err))
		}
	}
	RESULT solve(sizeType s,Vec& lambda,Eigen::Matrix<T,-1,-1>& U)
	{
		if(!_B)
		{
			setKrylovB(boost::shared_ptr<KrylovMatrix<T,KERNEL_TYPE> >(new IdentityKrylovMatrix<T,KERNEL_TYPE>(_A->n())));
		}
		checkMatrix();
		sizeType n=_A->n();
		//adjust 
		{
			s=std::min<sizeType>(s,_A->n());
			INFOV("User want %d eigenpairs",s)
			lambda.resize(s);
			U.resize(n,s);
			_iterCount.resize(s);
		}
		//compute
		for(sizeType eig=1;eig<=s;eig++)
		{
			INFOV("Computing the %dth pair",eig)
			findInit(U,eig);
			if(_cb)_cb->reset();

			T qk,eta;
			sizeType iter=0;
			while(true)
			{
				RESULT res=solve(U,eig,iter,qk,eta);
				if(res == Solver<T,KERNEL_TYPE>::SUCCESSFUL)
					break;
				else if(res == Solver<T,KERNEL_TYPE>::NOT_CONVERGENT)
					return res;
				_xk/=KERNEL_TYPE::norm(_xk);
			}

			_iterCount[eig-1]=iter;
			lambda[eig-1]=qk;
			U.col(eig-1)=_xk/std::sqrt(eta);
		}
		return Solver<T,KERNEL_TYPE>::SUCCESSFUL;
	}
	RESULT solve(const Eigen::Matrix<T,-1,-1>& U,sizeType eig,sizeType& iter,T& qk,T& eta)
	{
		T betak,gamma,qkLast=0.0f,alphak;
		T gk_gkM_Last,gk_gkM;
		T thres1,thres2;
		//step 1
		betak=0.0f;
		KERNEL_TYPE::copy(_xk,_pk);
		projectOut(U,eig,_pk,_xk);
		KERNEL_TYPE::zero(_pk);
		//step 2
		_A->multiply(_xk,_xkA);
		_B->multiply(_xk,_xkB);
		gamma=KERNEL_TYPE::dot(_xk,_xkA);
		eta=KERNEL_TYPE::dot(_xk,_xkB);
		qk=gamma/eta;
		//step
		sizeType innerIter=0;
		for(;innerIter<_maxIterInner && iter<_maxIter;innerIter++,iter++)
		{
			//3.1
			KERNEL_TYPE::copy(_xkA,_gk);
			KERNEL_TYPE::addScaled(-qk,_xkB,_gk);
			{
				//3.11 part2
				thres2=KERNEL_TYPE::norm(_gk)/std::sqrt(gamma);
				if(_cb)_cb->operator()(_xk,thres2,iter);
				if(thres2 < _eps2)
					break;
			}
			KERNEL_TYPE::scale(2.0f/eta,_gk);
			//3.2
			T gk_gkMLast=KERNEL_TYPE::dot(_gk,_gkM);
			_pre.solve(_gk,_gkM);
			//3.3
			gk_gkM=KERNEL_TYPE::dot(_gk,_gkM);
			if(innerIter>0)
				betak=(gk_gkM-gk_gkMLast)/gk_gkM_Last;
			gk_gkM_Last=gk_gkM;
			//3.4
			KERNEL_TYPE::copy(_gkM,_pkTilde);
			KERNEL_TYPE::addScaled(betak,_pk,_pkTilde);
			//3.5
			projectOut(U,eig,_pkTilde,_pk);
			//3.6
			_A->multiply(_pk,_pkA);
			_B->multiply(_pk,_pkB);
			//3.7
			{
				T xAp=KERNEL_TYPE::dot(_pk,_xkA);
				T pAp=KERNEL_TYPE::dot(_pk,_pkA);
				T xBp=KERNEL_TYPE::dot(_pk,_xkB);
				T pBp=KERNEL_TYPE::dot(_pk,_pkB);

				T a=xBp*pAp-xAp*pBp;
				T b=eta*pAp-gamma*pBp;
				T c=eta*xAp-gamma*xBp;
				T maxDenom=std::max(std::abs(a),std::max(std::abs(b),std::abs(c)));
				a/=maxDenom;
				b/=maxDenom;
				c/=maxDenom;

				T delta=b*b-4.0f*a*c;
				if(b > 0.0f)
					alphak=-2.0f*c/(b+sqrt(delta));
				else alphak=(-b+sqrt(delta))/(2.0f*a);
				
				//efficient method for updating gamma and eta,
				//but we use more stable version below
#define CHEAP_UPDATE_GAMMA_ETA
#ifdef CHEAP_UPDATE_GAMMA_ETA
				gamma=gamma+2.0f*xAp*alphak+pAp*alphak*alphak;
				eta=eta+2.0f*xBp*alphak+pBp*alphak*alphak;
#endif
			}
			//3.8
			KERNEL_TYPE::addScaled(alphak,_pk,_xk);
			KERNEL_TYPE::addScaled(alphak,_pkA,_xkA);
			KERNEL_TYPE::addScaled(alphak,_pkB,_xkB);
#ifndef CHEAP_UPDATE_GAMMA_ETA
			gamma=KERNEL_TYPE::dot(_xkA,_xk);
			eta=KERNEL_TYPE::dot(_xkB,_xk);
#endif
			//3.9
			if(gamma > 1E30f || eta > 1E30f)
				return Solver<T,KERNEL_TYPE>::USER_REQUEST;
			//3.10
			qk=gamma/eta;
			//INFOV("Update qk: %f",qk)
			if(qk < 0.0f)
			{
				if(gamma < 0.0f)
					NOTIFY_MSG("A not positive definite!")
				if(eta < 0.0f)
					NOTIFY_MSG("B not positive definite!")
				exit(EXIT_FAILURE);
			}
			//3.11 part1
			thres1=std::abs(qk-qkLast)/qk;
			if(thres1 < _eps1)
				return Solver<T,KERNEL_TYPE>::USER_REQUEST;//return Solver<T,KERNEL_TYPE>::NOT_CONVERGENT;
			qkLast=qk;
		}
		if(iter == _maxIter)
			return Solver<T,KERNEL_TYPE>::NOT_CONVERGENT;
		if(innerIter == _maxIterInner)
			return Solver<T,KERNEL_TYPE>::USER_REQUEST;
		return Solver<T,KERNEL_TYPE>::SUCCESSFUL;
	}
	void findInit(const Eigen::Matrix<T,-1,-1>& U,sizeType nr)
	{
		Vec tmp;
		KERNEL_TYPE::copy(_xk,tmp);
		while(true)
		{
			tmp.setRandom();
			projectOut(U,nr,tmp,_xk);
			if(KERNEL_TYPE::norm(_xk) > 1E-3f)
				break;
		}
		//check orthogonality
		_B->multiply(_xk,tmp);
		for(sizeType i=0;i<nr-1;i++)
			INFOV("\tOrthogonal residual: %f",U.col(i).dot(tmp))
	}
	void projectOut(const Eigen::Matrix<T,-1,-1>& U,sizeType nr,const Vec& from,Vec& to,bool includeU0=true) const
	{
		Vec tmp;
		KERNEL_TYPE::copy(from,tmp);
		KERNEL_TYPE::copy(from,to);

		_B->multiply(from,tmp);
		if(includeU0)
		for(sizeType i=0;i<_U0.cols();i++)
			to-=_U0.col(i).dot(tmp)*_U0.col(i);
		for(sizeType i=0;i<nr-1;i++)
			to-=U.col(i).dot(tmp)*U.col(i);
	}
	const vector<sizeType>& getIterationsCount() const{return _iterCount;}
	void checkMatrix() const
	{
		if(!_A || !_B)
		{
			NOTIFY_MSG("Matrix A or B not provieded!")
			exit(EXIT_FAILURE);
		}
		if(_A->n() != _B->n())
		{
			NOTIFY_MSG("Matrix A and B are of different size!")
			exit(EXIT_FAILURE);
		}
	}
protected:
	PRECONDITIONER_TYPE _pre;
	boost::shared_ptr<KrylovMatrix<T,KERNEL_TYPE> > _A,_B;
	Eigen::Matrix<T,-1,-1> _U0;
    boost::shared_ptr<Callback<T,KERNEL_TYPE> > _cb;
	//data
	T _eps1,_eps2;
	sizeType _maxIter,_maxIterInner;
	//temporary
	Vec _xk,_xkA,_xkB;
	Vec _pk,_pkA,_pkB,_pkTilde;
	Vec _gk,_gkM;
	vector<sizeType> _iterCount;
};

PRJ_END

#endif