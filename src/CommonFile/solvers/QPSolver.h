#ifndef QP_SOLVER_H
#define QP_SOLVER_H

#include "solvers/LinearSolver.h"

PRJ_BEGIN

template <typename T,typename KERNEL_TYPE=Kernel<T> >
class MPRGPQPSolver;
//a set of trivial preconditioner
template <typename T,typename KERNEL_TYPE=Kernel<T> >
class InFaceNoPreconSolver
{
public:
	typedef typename KERNEL_TYPE::Vec Vec;
	virtual void solve(const Vec& rhs,Vec& result,const std::vector<char>& face)
	{
		MPRGPQPSolver<T,KERNEL_TYPE>::MASK_FACE(rhs,result,face);
	}
};
template <typename T,typename KERNEL_TYPE=Kernel<T> >
class DiagonalInFacePreconSolver : public InFaceNoPreconSolver<T,KERNEL_TYPE>
{
public:
    using typename InFaceNoPreconSolver<T,KERNEL_TYPE>::Vec;
	DiagonalInFacePreconSolver(sizeType n):_invDiag(n){}
	DiagonalInFacePreconSolver(const FixedSparseMatrix<T,KERNEL_TYPE>& m)
	{
		_invDiag.resize(m.rows());
		for(sizeType i=0;i<m.rows();i++)
			_invDiag[i]=1.0f/m(i,i);
	}
	DiagonalInFacePreconSolver(const Eigen::Matrix<T,-1,-1>& m)
	{
		_invDiag=m.diagonal();
		_invDiag.array()=1.0f/_invDiag.array();
	}
	virtual void solve(const Vec& rhs,Vec& result,const std::vector<char>& face)
	{
		result.array()=rhs.array()*_invDiag.array();
		MPRGPQPSolver<T,KERNEL_TYPE>::MASK_FACE(rhs,result,face);
	}
	Vec _invDiag;
};
//use Preconditioned MPRGP as outter iteration (R-linear, large scale, simple box constraint only)
template <typename T,typename KERNEL_TYPE>
class MPRGPQPSolver : public Solver<T,KERNEL_TYPE>
{
public:
	typedef typename KERNEL_TYPE::Vec Vec;
	//constructor
	virtual ~MPRGPQPSolver(){}
	void reset(const Vec& L,const Vec& H)
	{
		_L=&L;
		_H=&H;
		setSolverParameters(1e-5f,1000);
		KERNEL_TYPE::copy(L,_g);
		KERNEL_TYPE::copy(L,_p);
		KERNEL_TYPE::copy(L,_z);
		KERNEL_TYPE::copy(L,_beta);
		KERNEL_TYPE::copy(L,_phi);
		KERNEL_TYPE::copy(L,_gp);
		_face.resize(L.rows());
		_pre.reset(new InFaceNoPreconSolver<T,KERNEL_TYPE>);
		_Gamma=1.0f;
	}
    virtual void setSolverParameters(T toleranceFactor,sizeType maxIterations) {
        _maxIterations=maxIterations;
		_toleranceFactor=toleranceFactor;
        if(_toleranceFactor<1e-30f)
            _toleranceFactor=1e-30f;
    }
	virtual void setMatrix(const FixedSparseMatrix<T,KERNEL_TYPE>& A,bool syncPrecon)
	{
		_A.reset(new DefaultKrylovMatrix<T>(A));
		_pre.reset(new DiagonalInFacePreconSolver<T>(A));
		_alphaBar=2.0f/specRad(*_A);
	}
	virtual void setMatrix(const Eigen::Matrix<T,-1,-1>& A)
	{
		_A.reset(new DefaultDenseEigenKrylovMatrix<T>(A));
		_pre.reset(new DiagonalInFacePreconSolver<T>(A));
		_alphaBar=2.0f/specRad(*_A);
	}
    virtual void setKrylovMatrix(boost::shared_ptr<KrylovMatrix<T,KERNEL_TYPE> > matrix)
	{
		_A=matrix;
		_alphaBar=2.0f/specRad(*_A);
	}
	virtual void setPre(boost::shared_ptr<InFaceNoPreconSolver<T,KERNEL_TYPE> > pre){_pre=pre;}
	//methods
	static T specRad(const KrylovMatrix<T>& G,Vec* ev=NULL,const T& eps=1E-3f)
	{
		T delta;
		Vec tmp,tmpOut;
		tmp.resize(G.n());
		tmpOut.resize(G.n());
		tmp.setRandom();
		tmp.normalize();

		//power method
		for(sizeType iter=0;;iter++)
		{
			G.multiply(tmp,tmpOut);
			T normTmpOut=tmpOut.norm();
			if(normTmpOut < ScalarUtil<T>::scalar_eps){
				if(ev)*ev=tmp;
				return ScalarUtil<T>::scalar_eps;
			}

			tmpOut/=normTmpOut;
			delta=(tmpOut-tmp).norm();
			//INFOV("Power Iter %d Err: %f, SpecRad: %f",iter,delta,normTmpOut)
			if(delta <= eps){
				if(ev)*ev=tmp;
				return normTmpOut;
			}
			tmp=tmpOut;
		}
	}
	const std::vector<char>& getFace() const{return _face;}
    bool checkKKT(const Vec& result,const Vec& B)
	{
		_A->multiply(result,_g);
        KERNEL_TYPE::sub(_g,B,_g);
		for(sizeType i=0;i<result.size();i++)
		{
            if(std::abs(result[i]-(*_L)[i]) < _toleranceFactor){
				if(_g[i] < -ScalarUtil<T>::scalar_eps)
					return false;
            }else if(std::abs(result[i]-(*_H)[i]) < _toleranceFactor){
				if(_g[i] > ScalarUtil<T>::scalar_eps)
					return false;
            }else if(std::abs(_g[i]) > _toleranceFactor)
				return false;
		}
		return true;
	}
    virtual typename Solver<T,KERNEL_TYPE>::SOLVER_RESULT solve(const Vec& rhs,Vec& result)
	{
		//declaration
		T alphaCG,alphaF,beta;
		Vec& AP=_gp;		//clever reuse ?
		Vec& y=_beta;		//clever reuse ?
		Vec& xTmp=_beta;	//clever reuse ?
		Vec& D=_phi;		//clever reuse ?

		//initialize
		_A->multiply(result,_g);
		KERNEL_TYPE::sub(_g,rhs,_g);
		DECIDE_FACE(result,*_L,*_H,_face);
		_pre->solve(_g,_z,_face);
		KERNEL_TYPE::copy(_z,_p);	//_p = phi(x)
        if(Solver<T,KERNEL_TYPE>::_cb)Solver<T,KERNEL_TYPE>::_cb->reset();

		//MPRGP iteration
		sizeType iteration;
		for(iteration=0;iteration<_maxIterations;iteration++)
		{
			//test termination
			MASK_FACE(_g,_phi,_face);
			BETA(_g,_beta,_face);
			KERNEL_TYPE::add(_phi,_beta,_gp);
            Solver<T,KERNEL_TYPE>::_residualOut=KERNEL_TYPE::norm(_gp);
            if(Solver<T,KERNEL_TYPE>::_cb)
                (*Solver<T,KERNEL_TYPE>::_cb)(result,Solver<T,KERNEL_TYPE>::_residualOut,iteration);
            if(Solver<T,KERNEL_TYPE>::_residualOut < _toleranceFactor){
                Solver<T,KERNEL_TYPE>::_iterationsOut=iteration;
				return Solver<T,KERNEL_TYPE>::SUCCESSFUL;
			}

			//test proportional x
			//Note in the box constrained version of MPRGP
			//Test condition beta*beta <= gamma*gamma*phi*phiTilde
			//Must be replaced with beta*betaTilde <= gamma*gamma*phi*phiTilde
			//This is because proportioning step can also be blocked 
			//So we also introduce proportioning expansion step, see below
			if(BETATBETA(result,_g,*_L,*_H,_alphaBar,_beta,_face) <= _Gamma*_Gamma*PHITPHI(result,*_L,*_H,_alphaBar,_phi)){
				
				//prepare conjugate gradient
				_A->multiply(_p,AP);
				alphaCG=KERNEL_TYPE::dot(_g,_z)/
						KERNEL_TYPE::dot(AP,_p);
				KERNEL_TYPE::copy(result,y);
				KERNEL_TYPE::addScaled(-alphaCG,_p,y);
				alphaF=stepLimit(result,_p);

				if(alphaCG < alphaF){
					//conjugate gradient step
                    if(Solver<T,KERNEL_TYPE>::_cb)INFO("\tConjugate Gradient Step")
					KERNEL_TYPE::copy(y,result);
					KERNEL_TYPE::addScaled(-alphaCG,AP,_g);
					_pre->solve(_g,_z,_face);	//face must not change, in conjugate gradient step
					beta=KERNEL_TYPE::dot(AP,_z)/
						 KERNEL_TYPE::dot(AP,_p);
					KERNEL_TYPE::scale(-beta,_p);
					KERNEL_TYPE::add(_z,_p,_p);
				}else{
					//expansion step
                    if(Solver<T,KERNEL_TYPE>::_cb)INFO("\tExpansion Step")
					KERNEL_TYPE::copy(result,xTmp);
					KERNEL_TYPE::addScaled(-alphaF,_p,xTmp);
					KERNEL_TYPE::addScaled(-alphaF,AP,_g);
					DECIDE_FACE(xTmp,*_L,*_H,_face);MASK_FACE(_g,_phi,_face);	//decide face for xTmp
					KERNEL_TYPE::addScaled(-_alphaBar,_phi,xTmp);
					project(xTmp,result);
					//restart CG
					_A->multiply(result,_g);
					KERNEL_TYPE::sub(_g,rhs,_g);
					DECIDE_FACE(result,*_L,*_H,_face);
					_pre->solve(_g,_z,_face);
					KERNEL_TYPE::copy(_z,_p);	//face can change
				}

			}else{

				//prepare proportioning
				KERNEL_TYPE::copy(_beta,D);	//not needed for clever reused version
				_A->multiply(D,AP);
				alphaCG=KERNEL_TYPE::dot(_g,D)/
						KERNEL_TYPE::dot(AP,D);
				alphaF=stepLimit(result,D);
				
				if(alphaCG < alphaF){
					//proportioning step
                    if(Solver<T,KERNEL_TYPE>::_cb)INFO("\tProportioning Step")
					KERNEL_TYPE::addScaled(-alphaCG,D,result);
					//restart CG
					KERNEL_TYPE::addScaled(-alphaCG,AP,_g);
					DECIDE_FACE(result,*_L,*_H,_face);
					_pre->solve(_g,_z,_face);
					KERNEL_TYPE::copy(_z,_p);	//face can change
				}else{
					//proportioning expansion step
                    if(Solver<T,KERNEL_TYPE>::_cb)INFO("\tProportioning Expansion Step")
					KERNEL_TYPE::copy(result,xTmp);
					KERNEL_TYPE::addScaled(-alphaF,D,xTmp);
					KERNEL_TYPE::addScaled(-alphaF,AP,_g);
					DECIDE_FACE(xTmp,*_L,*_H,_face);BETA(_g,D,_face);	//reused D
					KERNEL_TYPE::addScaled(-_alphaBar,D,xTmp);
					project(xTmp,result);
					//restart CG
					_A->multiply(result,_g);
					KERNEL_TYPE::sub(_g,rhs,_g);
					DECIDE_FACE(result,*_L,*_H,_face);
					_pre->solve(_g,_z,_face);
					KERNEL_TYPE::copy(_z,_p);	//face can change
				}

			}
		}

        Solver<T,KERNEL_TYPE>::_iterationsOut=iteration;
        return Solver<T,KERNEL_TYPE>::NOT_CONVERGENT;
	}
	static void MASK_FACE(const Vec& in,Vec& out,const std::vector<char>& face)
	{
		OMP_PARALLEL_FOR_
		for(sizeType i=0;i<in.rows();i++)
			if(face[i] != 0)
				out[i]=0.0f;
			else out[i]=in[i];
	}
protected:
	T stepLimit(const Vec& X,const Vec& D) const
	{
		T ret=ScalarUtil<T>::scalar_max;
		T tmp;
		#pragma omp parallel private(tmp)
		{
			tmp=ScalarUtil<T>::scalar_max;
			#pragma omp for
			for(sizeType i=0;i<_A->n();i++)
			{
				if(D[i] > ScalarUtil<T>::scalar_eps && X[i] > (*_L)[i])	//handle rounding err
					tmp=std::min<T>(tmp,(X[i]-(*_L)[i])/D[i]);
				else if(D[i] < -ScalarUtil<T>::scalar_eps && X[i] < (*_H)[i])	//handle rounding err
					tmp=std::min<T>(tmp,(X[i]-(*_H)[i])/D[i]);
			}

			OMP_CRITICAL_
			ret=std::min<T>(ret,tmp);
		}
		return ret;
	}
	void project(const Vec& in,Vec& out) const
	{
		OMP_PARALLEL_FOR_
		for(sizeType i=0;i<_A->n();i++)
			out[i]=std::min<T>(std::max<T>(in[i],(*_L)[i]),(*_H)[i]);
	}
	static void BETA(const Vec& in,Vec& out,const std::vector<char>& face)
	{
		OMP_PARALLEL_FOR_
		for(sizeType i=0;i<in.rows();i++)
			if(face[i] == 0)
				out[i]=0.0f;
			else if(face[i] == 1)
				out[i]=std::max<T>(in[i],0.0f);
			else out[i]=std::min<T>(in[i],0.0f);
	}
	static void DECIDE_FACE(const Vec& x,const Vec& L,const Vec& H,std::vector<char>& face)
	{
		face.assign(x.rows(),0);
		OMP_PARALLEL_FOR_
		for(sizeType i=0;i<x.rows();i++)
            if(std::abs(x[i]-L[i]) < ScalarUtil<T>::scalar_eps)
				face[i]=-1;
            else if(std::abs(x[i]-H[i]) < ScalarUtil<T>::scalar_eps)
				face[i]=1;
	}
	static T PHITPHI(const Vec& x,const Vec& L,const Vec& H,const T& alphaBar,const Vec& phi)
	{
		T phiTphi=0.0f;
		#ifdef _MSC_VER
        OMP_PARALLEL_FOR_I(OMP_ADD(phiTphi))
		#endif
		for(sizeType i=0;i<x.rows();i++){
			T phiTilde=0.0f;
			if(phi[i] > 0.0f && x[i] > L[i])	//handle rounding error
				phiTilde=std::min<T>((x[i]-L[i])/alphaBar,phi[i]);
			else if(phi[i] < 0.0f && x[i] < H[i])	//handle rounding error
				phiTilde=std::max<T>((x[i]-H[i])/alphaBar,phi[i]);

			ASSERT(phiTilde*phi[i] >= 0.0f)
			phiTphi+=phiTilde*phi[i];
		}
		return phiTphi;
	}
	static T BETATBETA(const Vec& x,const Vec& g,const Vec& L,const Vec& H,const T& alphaBar,const Vec& beta,std::vector<char>& face)
	{
		T betaTbeta=0.0f;
		#ifdef _MSC_VER
        OMP_PARALLEL_FOR_I(OMP_ADD(betaTbeta))
		#endif
		for(sizeType i=0;i<x.rows();i++)
		{
			T betaTilde=0.0f;
			if(face[i] == -1 && g[i] < 0.0f && x[i] < H[i])	//handle rounding error
				betaTilde=std::max<T>((x[i]-H[i])/alphaBar,g[i]);
			else if(face[i] == 1 && g[i] > 0.0f && x[i] > L[i])	//handle rounding error
				betaTilde=std::min<T>((x[i]-L[i])/alphaBar,g[i]);
			
			ASSERT(betaTilde*beta[i] >= 0.0f)
			betaTbeta+=betaTilde*beta[i];
		}
		return betaTbeta;
	}
protected:
	// internal structures
    boost::shared_ptr<InFaceNoPreconSolver<T,KERNEL_TYPE> > _pre;
	//problem
	boost::shared_ptr<KrylovMatrix<T> > _A;
	const Vec* _L;
	const Vec* _H;
	//parameter
	sizeType _maxIterations;
    T _toleranceFactor;
	T _Gamma,_alphaBar;
	//temporary
	std::vector<char> _face;
	Vec _g,_p,_z,_beta,_phi,_gp;
};

PRJ_END

#endif
