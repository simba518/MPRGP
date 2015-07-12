#ifndef _QPSOLVER_H_
#define _QPSOLVER_H_

#include <stdio.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <boost/shared_ptr.hpp>
#include <vector>

#include <omp.h>
#ifdef _MSC_VER
#define STRINGIFY(X) X
#define PRAGMA __pragma
#else
#define STRINGIFY(X) #X
#define PRAGMA _Pragma
#endif

#define OMP_PARALLEL_FOR_ PRAGMA(STRINGIFY(omp parallel for num_threads(OmpSettings::getOmpSettings().nrThreads()) schedule(dynamic,OmpSettings::getOmpSettings().szChunk())))
#define OMP_CRITICAL_ PRAGMA(STRINGIFY(omp critical))

typedef int sizeType;

template <typename T>
struct ScalarUtil;
template <>
struct ScalarUtil<float> {
  static float scalar_max;
  static float scalar_eps;
};
template <>
struct ScalarUtil<double> {
  static double scalar_max;
  static double scalar_eps;
};

template <typename T>
struct FixedSparseMatrix:public Eigen::SparseMatrix<T>{
public:
  template <typename VEC,typename VEC_OUT>
  void multiply(const VEC& x,VEC_OUT& result)const{
	result = (*this)*x;
  }
  double diag(const int i)const{
	return (*this)(i,i);
  }
  double funValue(const VectorXd &x)const{
	printf("unimplemented function is called\n");
	return 0.0f;
  }
};

// Iteration callback
template <typename T, typename FUNCTION = FixedSparseMatrix<T> >
struct Callback {
public:
  typedef Eigen::Matrix<T,-1,1> Vec;
  Callback(FUNCTION *func=NULL):function(func){}
  virtual ~Callback() {}
  void setFunction(FUNCTION &func){
	this->function = &func;
  }
  virtual void reset() {}
  virtual sizeType operator()(const Vec& x,const T& residueNorm,const sizeType& n) {
	printf("ITER %d: r=%f",(int)n,residueNorm);
	if ( function && 0 == n%10){
	  printf(", fun=%f", function->funValue(x));
	}
	return 0;
  }
  virtual sizeType operator()(const Vec& x,const Vec& g,const T& fx,const T& xnorm,const T& gnorm,const T& step,const sizeType& k,const sizeType& ls) {
	printf("ITER %d: f=%f, x=%f, g=%f, s=%f",(int)k,fx,xnorm,gnorm,step);
	return 0;
  }
  
private:
  FUNCTION *function;
};

template <typename T, typename FUNCTION>
struct Solver {
public:
  enum SOLVER_RESULT {
	SUCCESSFUL,
	NOT_CONVERGENT,
	USER_REQUEST,
	NEGATIVE_EIGENVALUE,
  };
  Solver():_cb((Callback<T, FUNCTION>*)NULL) {}
  virtual void setCallback(typename boost::shared_ptr<Callback<T, FUNCTION> >cb){_cb=cb;}

protected:
  sizeType _iterationsOut;
  T _residualOut;
  boost::shared_ptr<Callback<T, FUNCTION> > _cb;
};

template <typename T,typename MAT, typename KERNEL_TYPE, typename FUNCTION >
class MPRGPQPSolver;

//a set of trivial preconditioner
template <typename T,typename KERNEL_TYPE, typename MAT, typename FUNCTION>
class InFaceNoPreconSolver{

  typedef Eigen::Matrix<T,-1,1> Vec;

public:
  InFaceNoPreconSolver(const std::vector<char>& face):_face(face){}
  virtual typename Solver<T,FUNCTION>::SOLVER_RESULT solve(const Vec&rhs,Vec&result){
  	MPRGPQPSolver<T,MAT,KERNEL_TYPE, FUNCTION>::MASK_FACE(rhs,result,_face);
  	return Solver<T,FUNCTION>::SUCCESSFUL;
  }
  virtual void setMatrix(const MAT& matrix){}

protected:
  const std::vector<char>& _face;
};

template <typename T,typename KERNEL_TYPE, typename MAT, typename FUNCTION >
class DiagonalInFacePreconSolver : public InFaceNoPreconSolver<T,KERNEL_TYPE,MAT, FUNCTION>{

  typedef Eigen::Matrix<T,-1,1> Vec;

public:
  DiagonalInFacePreconSolver(const std::vector<char>& face)
	:InFaceNoPreconSolver<T,KERNEL_TYPE,MAT, FUNCTION>(face){}
  virtual typename Solver<T, FUNCTION>::SOLVER_RESULT solve(const Vec&rhs,Vec&result){
	for(sizeType i=0;i<rhs.rows();i++){
	  if(InFaceNoPreconSolver<T,KERNEL_TYPE,MAT,FUNCTION>::_face[i] == 0)
		result[i]=rhs[i]/_matrix->diag(i);
	  else result[i]=0.0f;
	}
	return Solver<T, FUNCTION>::SUCCESSFUL;
  }
  virtual void setMatrix(const MAT& matrix){
	_matrix=&matrix;
  }
  const MAT* _matrix;
};

// Vector Interface
template <typename T>
struct Kernel {
  template<typename T2>
  struct Rebind {
	typedef Kernel<T2> value;
  };
  typedef Eigen::Matrix<T,-1,1> Vec;
  static inline T dot(const Vec& x, const Vec& y,sizeType n=-1) {
	if(n == -1)return x.dot(y);
	else return x.block(0,0,n,1).dot(y.block(0,0,n,1));
  }
  static inline T norm(const Vec& x,sizeType n=-1) {
	if(n == -1)return x.norm();
	else return x.block(0,0,n,1).norm();
  }
  static inline void copy(const Vec& x,Vec& y,sizeType n=-1) {
	if(n == -1)y=x;
	else y.block(0,0,n,1)=x.block(0,0,n,1);
  }
  static inline void ncopy(const Vec& x,Vec& y,sizeType n=-1) {
	if(n == -1)y=-x;
	else y.block(0,0,n,1)=-x.block(0,0,n,1);
  }
  static inline sizeType indexAbsMax(const Vec& x,sizeType n=-1) {
	if(n == -1)n=x.size();
	sizeType maxInd = 0;
	T maxValue = (T)0.0f;
	for(sizeType i = 0; i < n; ++i) {
	  if(fabs(x[i]) > maxValue) {
		maxValue = fabs(x[i]);
		maxInd = i;
	  }
	}
	return maxInd;
  }
  static inline T absMax(const Vec& x,sizeType n=-1) {
	if(n == -1)return std::fabs(x[indexAbsMax(x)]);
	else return std::fabs(x.block(0,0,n,1)[indexAbsMax(x.block(0,0,n,1))]);
  }
  static inline void scale(T alpha,Vec& y,sizeType n=-1) {
	if(n == -1)y*=alpha;
	else y.block(0,0,n,1)*=alpha;
  }
  static inline void addScaled(const T& alpha, const Vec& x, Vec& y,sizeType n=-1) {
	if(n == -1)y+=x*alpha;
	else y.block(0,0,n,1)+=x.block(0,0,n,1)*alpha;
  }
  static inline void add(const Vec& x,const Vec& y,Vec& result,sizeType n=-1) {
	if(n == -1)result=x+y;
	else result.block(0,0,n,1)=x.block(0,0,n,1)+y.block(0,0,n,1);
  }
  static inline void sub(const Vec& x,const Vec& y,Vec& result,sizeType n=-1) {
	if(n == -1)result=x-y;
	else result.block(0,0,n,1)=x.block(0,0,n,1)-y.block(0,0,n,1);
  }
  static inline void zero(Vec& v,sizeType n=-1) {set(v,0.0f,n);}
  static inline void set(Vec& v,const T& val,sizeType n=-1) {
	if(n == -1)v.setConstant(val);
	else v.block(0,0,n,1).setConstant(val);
  }
  template <typename A,typename B>
  static inline void addMulT(A& x,const T& b,const B& c){x+=b*c;}
  template <typename A,typename B>
  static inline void subMulT(A& x,const T& b,const B& c){x-=b*c;}
};

// solve for 
// 
// min 1/2*x^t*A*x-x^t*B, s.t.  L<=x<=H.
// 
// use Preconditioned MPRGP as outter iteration (R-linear, large scale, simple
// box constraint only)
template <typename T,typename MAT=FixedSparseMatrix<T>, typename KERNEL_TYPE=Kernel<T>, typename FUNCTION = FixedSparseMatrix<T> >
class MPRGPQPSolver:public Solver<T, FUNCTION>{

public:
  typedef Eigen::Matrix<T,-1,1> Vec;
  // constructor
  MPRGPQPSolver(const MAT& A,const Vec& B,const Vec& L,const Vec& H,const bool precond=true)
  	:_A(A),_B(B),_L(L),_H(H){
  	setSolverParameters(1e-5f,1000);
  	KERNEL_TYPE::copy(_B,_g);
  	KERNEL_TYPE::copy(_B,_p);
  	KERNEL_TYPE::copy(_B,_z);
  	KERNEL_TYPE::copy(_B,_beta);
  	KERNEL_TYPE::copy(_B,_phi);
  	KERNEL_TYPE::copy(_B,_gp);
  	_face.resize(_B.rows());
	if (precond){
	  _pre.reset(new DiagonalInFacePreconSolver<T,KERNEL_TYPE,MAT,FUNCTION>(_face));
	  _pre->setMatrix(_A);
	}else{
	  _pre.reset(new InFaceNoPreconSolver<T,KERNEL_TYPE,MAT,FUNCTION>(_face));
	}
  }
  virtual ~MPRGPQPSolver(){}
  virtual void setSolverParameters(T toleranceFactor,sizeType maxIterations) {
	_maxIterations=maxIterations;
	_toleranceFactor=toleranceFactor;
	if(_toleranceFactor<1e-30f)
	  _toleranceFactor=1e-30f;
	_Gamma=1.0f;
	_alphaBar=2.0f/specRad(_A);
  }
  virtual void setInFacePreconditioner(boost::shared_ptr<InFaceNoPreconSolver<T,KERNEL_TYPE,MAT,FUNCTION> > pre){
  	_pre=pre;
  	_pre->setMatrix(_A);
  }

  //methods
  static T specRad(const MAT& G,Vec* ev=NULL,const T& eps=1E-3f){

  	T delta;
  	Vec tmp,tmpOut;
  	tmp.resize(G.rows());
  	tmpOut.resize(G.rows());
  	tmp.setRandom();
  	tmp.normalize();

  	//power method
  	// for(sizeType iter=0;;iter++){ /// @todo
  	for(sizeType iter=0;iter <= 1000;iter++){
  	  G.multiply(tmp,tmpOut);
  	  T normTmpOut=tmpOut.norm();
  	  if(normTmpOut < ScalarUtil<T>::scalar_eps){
  		if(ev)*ev=tmp;
  		return ScalarUtil<T>::scalar_eps;
  	  }
  	  tmpOut/=normTmpOut;
  	  delta=(tmpOut-tmp).norm();
  	  printf("Power Iter %d Err: %f, SpecRad: %f\n",iter,delta,normTmpOut);
	  if(delta <= eps){
		if(ev)*ev=tmp;
		return normTmpOut;
	  }
  	  tmp=tmpOut;
  	}
  }
  const std::vector<char>& getFace() const{return _face;}
  bool checkKKT(const Vec& result){
  	_A.multiply(result,_g);
  	KERNEL_TYPE::sub(_g,_B,_g);
  	for(sizeType i=0;i<result.size();i++){
	  if(abs(result[i]-_L[i]) < ScalarUtil<T>::scalar_eps){
		if(_g[i] < -ScalarUtil<T>::scalar_eps)
		  return false;
	  }else if(abs(result[i]-_H[i]) < ScalarUtil<T>::scalar_eps){
		if(_g[i] > ScalarUtil<T>::scalar_eps)
		  return false;
	  }else if(abs(_g[i]) > ScalarUtil<T>::scalar_eps)
		return false;
	}
  	return true;
  }
  ///@note the initial value in result should be in the bound.
  virtual typename Solver<T,FUNCTION >::SOLVER_RESULT solve(Vec&result){

  	//declaration
  	T alphaCG,alphaF,beta;
  	Vec& AP=_gp;		//clever reuse ?
  	Vec& y=_beta;		//clever reuse ?
  	Vec& xTmp=_beta;	//clever reuse ?
  	Vec& D=_phi;		//clever reuse ?

  	//initialize
  	_A.multiply(result,_g);
  	KERNEL_TYPE::sub(_g,_B,_g);
  	DECIDE_FACE(result,_L,_H,_face);
  	_pre->solve(_g,_z);
  	KERNEL_TYPE::copy(_z,_p);	//_p = phi(x)
  	if(_cb)_cb->reset(); 

  	//MPRGP iteration
  	sizeType iteration;
  	for(iteration=0;iteration<_maxIterations;iteration++){

	  //test termination
	  MASK_FACE(_g,_phi,_face);
	  BETA(_g,_beta,_face);
	  KERNEL_TYPE::add(_phi,_beta,_gp);
	  _residualOut=KERNEL_TYPE::norm(_gp);
	  if(_cb) (*_cb)(result,_residualOut,iteration);
	  if(_residualOut < _toleranceFactor){
		_iterationsOut=iteration;
		return Solver<T, FUNCTION>::SUCCESSFUL;
	  }

	  //test proportional x
	  //Note in the box constrained version of MPRGP
	  //Test condition beta*beta <= gamma*gamma*phi*phiTilde
	  //Must be replaced with beta*betaTilde <= gamma*gamma*phi*phiTilde
	  //This is because proportioning step can also be blocked 
	  //So we also introduce proportioning expansion step, see below
	  if(BETATBETA(result,_g,_L,_H,_alphaBar,_beta,_face) <= _Gamma*_Gamma*PHITPHI(result,_L,_H,_alphaBar,_phi)){
				
		//prepare conjugate gradient
		_A.multiply(_p,AP);
		alphaCG=KERNEL_TYPE::dot(_g,_z)/
		  KERNEL_TYPE::dot(AP,_p);
		KERNEL_TYPE::copy(result,y);
		KERNEL_TYPE::addScaled(-alphaCG,_p,y);
		alphaF=stepLimit(result,_p);

		if(alphaCG < alphaF){
		  //conjugate gradient step
		  if(_cb) printf("\tConjugate Gradient Step\n");
		  KERNEL_TYPE::copy(y,result);
		  KERNEL_TYPE::addScaled(-alphaCG,AP,_g);
		  _pre->solve(_g,_z);	//face must not change, in conjugate gradient step
		  beta=KERNEL_TYPE::dot(AP,_z)/
			KERNEL_TYPE::dot(AP,_p);
		  KERNEL_TYPE::scale(-beta,_p);
		  KERNEL_TYPE::add(_z,_p,_p);
		}else{
		  //expansion step
		  if(_cb) printf("\tExpansion Step\n");
		  KERNEL_TYPE::copy(result,xTmp);
		  KERNEL_TYPE::addScaled(-alphaF,_p,xTmp);
		  KERNEL_TYPE::addScaled(-alphaF,AP,_g);
		  DECIDE_FACE(xTmp,_L,_H,_face);MASK_FACE(_g,_phi,_face);	//decide face for xTmp
		  KERNEL_TYPE::addScaled(-_alphaBar,_phi,xTmp);
		  project(xTmp,result);
		  //restart CG
		  _A.multiply(result,_g);
		  KERNEL_TYPE::sub(_g,_B,_g);
		  DECIDE_FACE(result,_L,_H,_face);
		  _pre->solve(_g,_z);
		  KERNEL_TYPE::copy(_z,_p);	//face can change
		}

	  }else{

		//prepare proportioning
		KERNEL_TYPE::copy(_beta,D);	//not needed for clever reused version
		_A.multiply(D,AP);
		alphaCG=KERNEL_TYPE::dot(_g,D)/
		  KERNEL_TYPE::dot(AP,D);
		alphaF=stepLimit(result,D);
				
		if(alphaCG < alphaF){
		  //proportioning step
		  if(_cb) printf("\tProportioning Step\n");
		  KERNEL_TYPE::addScaled(-alphaCG,D,result);
		  //restart CG
		  KERNEL_TYPE::addScaled(-alphaCG,AP,_g);
		  DECIDE_FACE(result,_L,_H,_face);
		  _pre->solve(_g,_z);
		  KERNEL_TYPE::copy(_z,_p);	//face can change
		}else{
		  //proportioning expansion step
		  if(_cb) printf("\tProportioning Expansion Step\n");
		  KERNEL_TYPE::copy(result,xTmp);
		  KERNEL_TYPE::addScaled(-alphaF,D,xTmp);
		  KERNEL_TYPE::addScaled(-alphaF,AP,_g);
		  DECIDE_FACE(xTmp,_L,_H,_face);BETA(_g,D,_face);	//reused D
		  KERNEL_TYPE::addScaled(-_alphaBar,D,xTmp);
		  project(xTmp,result);
		  //restart CG
		  _A.multiply(result,_g);
		  KERNEL_TYPE::sub(_g,_B,_g);
		  DECIDE_FACE(result,_L,_H,_face);
		  _pre->solve(_g,_z);
		  KERNEL_TYPE::copy(_z,_p);	//face can change
		}
	  }
	}

  	_iterationsOut=iteration;
  	return Solver<T, FUNCTION>::NOT_CONVERGENT;
  }
  static void MASK_FACE(const Vec& in,Vec& out,const std::vector<char>& face){
  	OMP_PARALLEL_FOR_
  	  for(sizeType i=0;i<in.rows();i++)
  		if(face[i] != 0)
  		  out[i]=0.0f;
  		else out[i]=in[i];
  }

protected:
  T stepLimit(const Vec& X,const Vec& D) const{

  	T ret=ScalarUtil<T>::scalar_max;
  	T tmp;
#pragma omp parallel private(tmp)
  	{
  	  tmp=ScalarUtil<T>::scalar_max;
#pragma omp for
  	  for(sizeType i=0;i<_A.rows();i++)
  		{
  		  if(D[i] > ScalarUtil<T>::scalar_eps && X[i] > _L[i])	//handle rounding err
  			tmp=std::min<T>(tmp,(X[i]-_L[i])/D[i]);
  		  else if(D[i] < -ScalarUtil<T>::scalar_eps && X[i] < _H[i])	//handle rounding err
  			tmp=std::min<T>(tmp,(X[i]-_H[i])/D[i]);
  		}

  	  OMP_CRITICAL_
  		ret=std::min<T>(ret,tmp);
  	}
  	return ret;
  }
  void project(const Vec& in,Vec& out) const{
  	OMP_PARALLEL_FOR_
  	  for(sizeType i=0;i<_A.rows();i++)
  		out[i]=std::min<T>(std::max<T>(in[i],_L[i]),_H[i]);
  }
  static void BETA(const Vec& in,Vec& out,const std::vector<char>& face){
  	OMP_PARALLEL_FOR_
  	  for(sizeType i=0;i<in.rows();i++)
  		if(face[i] == 0)
  		  out[i]=0.0f;
  		else if(face[i] == 1)
  		  out[i]=std::max<T>(in[i],0.0f);
  		else out[i]=std::min<T>(in[i],0.0f);
  }
  static void DECIDE_FACE(const Vec& x,const Vec& L,const Vec& H,std::vector<char>& face){
  	face.assign(x.rows(),0);
  	OMP_PARALLEL_FOR_
  	  for(sizeType i=0;i<x.rows();i++)
  		if(abs(x[i]-L[i]) < ScalarUtil<T>::scalar_eps)
  		  face[i]=-1;
  		else if(abs(x[i]-H[i]) < ScalarUtil<T>::scalar_eps)
  		  face[i]=1;
  }
  static T PHITPHI(const Vec& x,const Vec&L,const Vec&H,const T&alphaBar,const Vec&phi){

  	T phiTphi=0.0f;
#pragma omp parallel for reduction(+:phiTphi)
  	for(sizeType i=0;i<x.rows();i++){
  	  T phiTilde=0.0f;
  	  if(phi[i] > 0.0f && x[i] > L[i])	//handle rounding error
  		phiTilde=std::min<T>((x[i]-L[i])/alphaBar,phi[i]);
  	  else if(phi[i] < 0.0f && x[i] < H[i])	//handle rounding error
  		phiTilde=std::max<T>((x[i]-H[i])/alphaBar,phi[i]);
  	  assert(phiTilde*phi[i] >= 0.0f);
	  phiTphi+=phiTilde*phi[i];
  	}
  	return phiTphi;
  }
  static T BETATBETA(const Vec& x,const Vec& g,const Vec& L,const Vec& H,const T& alphaBar,const Vec& beta,std::vector<char>& face){

  	T betaTbeta=0.0f;
#pragma omp parallel for reduction(+:betaTbeta)
  	for(sizeType i=0;i<x.rows();i++){
	  T betaTilde=0.0f;
	  if(face[i] == -1 && g[i] < 0.0f && x[i] < H[i])	//handle rounding error
		betaTilde=std::max<T>((x[i]-H[i])/alphaBar,g[i]);
	  else if(face[i] == 1 && g[i] > 0.0f && x[i] > L[i])	//handle rounding error
		betaTilde=std::min<T>((x[i]-L[i])/alphaBar,g[i]);
	  assert(betaTilde*beta[i] >= 0.0f);
	  betaTbeta+=betaTilde*beta[i];
	}
  	return betaTbeta;
  }
protected:
  // internal structures
  boost::shared_ptr<InFaceNoPreconSolver<T, KERNEL_TYPE, MAT, FUNCTION> > _pre;
  //problem
  const MAT& _A;
  const Vec& _B;
  const Vec& _L;
  const Vec& _H;
  //parameter
  sizeType _maxIterations;
  T _toleranceFactor;
  T _Gamma,_alphaBar;
  //temporary
  std::vector<char> _face;
  Vec _g,_p,_z,_beta,_phi,_gp;

  using Solver<T,FUNCTION>::_residualOut;
  using Solver<T,FUNCTION>::_iterationsOut;
  using Solver<T,FUNCTION>::_cb;
};

#endif /* _QPSOLVER_H_ */
