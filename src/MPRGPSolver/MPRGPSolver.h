#ifndef _MPRGPSOLVER_H_
#define _MPRGPSOLVER_H_

#include <string>
#include "MPRGPPrecondition.h"
#include "MPRGPProjection.h"

namespace MATH{

  /**
   * frameworks for the MPRGP method, where the preconditioning and projecting 
   * method can be changed through the template parameters.
   * 
   */
  template <typename T, typename MAT, typename PROJECTOIN,typename PRECONDITION>
  class MPRGPBase{

	typedef Eigen::Matrix<T,-1,1> Vec;
	
  public:
	MPRGPBase(const MAT &A,const Vec &B,
			  PRECONDITION &precond, PROJECTOIN &projector,
			  const int max_it = 1000, const T tol=1e-3, Vec *ev = NULL):
	  A(A), B(B), precond(precond), projector(projector), name("MPRGP") {

	  setParameters(tol,max_it);
	  gamma=1.0f;
	  alpha_bar=2.0/specRad(A,ev,tol);
	  iterations_out = 0;
	  residual_out = 0.0f;
	  DEBUG_LOG(name << setprecision(16) << " alpha_bar: " << alpha_bar);
	}
	
	void setParameters(T toleranceFactor,size_t maxIterations) {
	  max_iterations=maxIterations;
	  tolerance_factor=toleranceFactor;
	  if(tolerance_factor<1e-30f)
		tolerance_factor=1e-30f;
	}
	size_t iterationsOut()const{return iterations_out;}
	T residualOut()const{return residual_out;}
	static T specRad(const MAT& G,Vec* ev=NULL,const T& eps=1E-3f){

	  T delta;
	  Vec tmp,tmpOut;
	  tmpOut.resize(G.rows());
	  if (ev != NULL && ev->size() == G.rows()){
		tmp = (*ev);
	  }else{
		tmp.resize(G.rows());
		tmp.setRandom(); ///@todo warm start ?
		tmp.normalize();
	  }

	  T normTmpOut = 1.0;
	  //power method
	  // for(size_t iter=0;;iter++){ /// @todo
	  size_t iter=0;
	  for(;iter <= 1000;iter++){
		G.multiply(tmp,tmpOut);
		normTmpOut=tmpOut.norm();
		if(normTmpOut < ScalarUtil<T>::scalar_eps){
		  if(ev)*ev=tmp;
		  normTmpOut = ScalarUtil<T>::scalar_eps;
		  break;
		}
		tmpOut/=normTmpOut;
		delta=(tmpOut-tmp).norm();
		// DEBUG_LOG(setprecision(10)<<"power delta: "<<delta);
		// printf("Power Iter %d Err: %f, SpecRad: %f\n",iter,delta,normTmpOut);
		if(delta <= eps){
		  if(ev)*ev=tmp;
		  break;
		}
		tmp=tmpOut;
	  }
	  INFO_LOG("power iter: "<<iter);
	  return normTmpOut;
	}

	void setName(const std::string name){
	  this->name = name;
	}

	bool writeVTK(const Vec &x, const std::string filename)const{

	  Vec points(x.size()*4);
	  points.head(x.size()) = x;
	  points.segment(x.size(), x.size()) = phi+x;
	  points.segment(x.size()*2, x.size()) = beta+x;
	  points.segment(x.size()*3, x.size()) = g+x;

	  ofstream out;
	  out.open(filename.c_str());
	  if (!out.is_open()){
		ERROR_LOG(name << " failed to open this file: "<<filename);
		return false;
	  }

	  // head
	  out << "# vtk DataFile Version 3.1 \n";	  
	  out << "write phi(x), beta(x) and g(x) \nASCII\nDATASET UNSTRUCTURED_GRID\n";

	  // points
	  out << "POINTS "<< points.size()/3  <<" FLOAT\n";
	  for (int i = 0; i < points.size(); i += 3){
		out << points[i+0] << " " << points[i+1] << " " << points[i+2] << "\n";
	  }
	  
	  // lines
	  const int xp = x.size()/3;
	  const int num_lines = 3*xp;
	  out << "CELLS " <<num_lines<< " "<< 3*num_lines << "\n";
	  for (int i = 0; i < xp; i++){
		out << 2 << " "<< i << " "<< i+xp << "\n";
		out << 2 << " "<< i << " "<< i+2*xp << "\n";
		out << 2 << " "<< i << " "<< i+3*xp << "\n";
	  }
	  
	  out << "CELL_TYPES "<<num_lines<<"\n";
	  for (int i = 0; i < num_lines; ++i){
		out << 3;
		if (i == num_lines-1) out << "\n"; else  out << " ";
	  }

	  const bool succ = out.good();
	  out.close();
	  return succ;
	}
	bool printFace()const{
	  const vector<char> &f = projector.getFace();
	  cout << "face: ";
	  for (int i = 0; i < f.size(); ++i)
		cout << (int)f[i] << " ";
	  cout << "\n";
	  return true;
	}
	void printSolveInfo()const{
	  INFO_LOG(name << " cg steps: "<< num_cg);
	  INFO_LOG(name << " exp steps: "<< num_exp);
	  INFO_LOG(name << " prop steps: "<< numprop);
	  INFO_LOG(name << " total iter: "<<iteration);
	  INFO_LOG(name << " residual: "<<residual_out);
	  INFO_LOG(name << " constraints: "<<COUNT_CONSTRAINTS(projector.getFace()));
	}

  protected:
	inline void initialize(const Vec &result){
	  
	  A.multiply(result,g);
	  g -= B;
	  projector.DECIDE_FACE(result);
	  projector.PHI(g, phi);
	  precond.solve(g,z, phi);     assert_eq(z,z);
	  assert_ge_ext(g.dot(z),0,"\ng.dot(phi) = "<<g.dot(phi)<<"\n|phi|="<<phi.norm()<<"\n|z|="<<z.norm());
	  p = z;
	  
	  result_code = -1;
	  num_cg = num_exp = numprop = iteration = 0;

	  debug_fun(fun_val = A.funcValue(result, B); cout << name << setprecision(12) << " func = " << fun_val<<endl;);
	}
	inline T computeGradients(const Vec &g, const bool comp_phi = true){

	  // debug_fun(assert(printFace()));
	  if (comp_phi)
		projector.PHI(g,phi);
	  projector.BETA(g,beta,phi);
	  gp = phi+beta;
	  const T residual = gp.norm();
	  // assert_le(phi.dot(beta),ScalarUtil<T>::scalar_eps*residual);
	  return residual;
	}
	inline bool proportional(const Vec &result)const{
	  
	  //test proportional x: beta*beta <= gamma*gamma*phi*phiTilde
	  const T left = beta.dot(beta);  assert_eq(left, left);
	  const T phitphi = projector.PHITPHI(result,alpha_bar,phi);
	  const T right = gamma*gamma*phitphi; assert_eq(right, right);
	  return left <= right;
	}
	inline void CGStep(const Vec &AP, T alpha_cg, Vec &result){

	  DEBUG_LOG(name << " cg_step");
	  num_cg ++;

	  assert_eq(alpha_cg, alpha_cg);
	  result -= alpha_cg*p;
	  g -= alpha_cg*AP;
	  projector.PHI(g, phi);  assert_ge(g.dot(phi),0);
	  precond.solve(g,z, phi);   assert_eq(z, z);  assert_ge(g.dot(z),0);
	  const T beta = (z.dot(AP)) / (p.dot(AP)); assert_eq(beta, beta);
	  p = z-beta*p;
	  debug_fun(fun_val = A.funcValue(result, B); cout << name << setprecision(12) << " func = " << fun_val<<endl;);
	}
	inline void ExpStep(const Vec &AP, T alpha_f, Vec &result){

	  DEBUG_LOG(name << " exp_step");
	  num_exp ++;

	  assert_eq(alpha_f, alpha_f);
	  Vec &xTmp=beta;
	  xTmp = result-alpha_f*p;
	  g -= alpha_f*AP;
	  projector.DECIDE_FACE(xTmp);
	  projector.PHI(g, phi);
	  xTmp -= alpha_bar*phi;
	  projector.project(xTmp,result);

	  A.multiply(result,g);
	  g -= B;
	  projector.DECIDE_FACE(result);
	  projector.PHI(g, phi); 
	  precond.solve(g,z, phi); assert_eq(z, z);
	  p = z;

	  debug_fun(fun_val = A.funcValue(result, B); cout << name << setprecision(12) << " func = " << fun_val<<endl;);
	}
	inline void PropStep(Vec &result){

	  DEBUG_LOG(name << " prop_step");
	  numprop ++;

	  const Vec &D = beta;
	  Vec& AD=gp;
	  A.multiply(D,AD);
	  const T ddad = D.dot(AD); 
	  assert_ne_ext(ddad,0,"|beta| = " << beta.norm()
					<< ", |phi| = " << phi.norm() 
					<< ", |face| = " << COUNT_CONSTRAINTS(projector.getFace()));
	  const T alpha_cg = (g.dot(D)) / ddad; assert_eq(alpha_cg, alpha_cg);
	  result -= alpha_cg*D;

	  Vec& xTmp=beta;
	  xTmp = result;
	  projector.project(xTmp,result);
	  g -= alpha_cg*AD;
	  projector.DECIDE_FACE(result);
	  projector.PHI(g, phi);
	  precond.solve(g,z, phi); assert_eq(z, z);
	  p = z;
	  debug_fun(fun_val = A.funcValue(result, B); cout<< name << setprecision(12) << " func = " << fun_val<<endl;);
	}

	inline T projectFunValue(const Vec &x)const{
	  static Vec tempx; /// @todo
	  projector.project(x, tempx);
	  return A.funcValue(tempx, B);
	}
	inline void ExpMonotonicStep(Vec &result){

	  DEBUG_LOG(name << " exp_step");
	  num_exp ++;

	  projector.project(result, y);
	  A.multiply(y,g);
	  g -= B;
	  projector.DECIDE_FACE(y);
	  projector.PHI(g, phi);
	  y -= alpha_bar*phi;
	  projector.project(y, result);

	  A.multiply(result,g);
	  g -= B;
	  projector.DECIDE_FACE(result);
	  projector.PHI(g, phi); 
	  precond.solve(g,z, phi); assert_eq(z, z);
	  p = z;

	  debug_fun(fun_val = A.funcValue(result, B); cout<< name << setprecision(12) << " func = " << fun_val<<endl;);
	}
	inline T CGMonotonicStep(Vec &result, bool &result_is_prop){

	  DEBUG_LOG(name << " begin cg_step");
	  assert( projector.isFeasible(result) );
	  Vec &AP=gp;
	  A.multiply(p, AP);
	  const T pd = p.dot(AP); assert_ge(pd, 0); // pd = p^t*A*p > 0
	  assert_ge_ext(g.dot(z),0,"\ng.dot(phi) = "<<g.dot(phi)<<"\n|phi|="<<phi.norm()<<"\n|beta|="<<beta.norm()<<"\n|z|="<<z.norm());
	  T alpha_cg = (z.dot(g)) / pd; assert_ge(alpha_cg, 0);

	  y = result - alpha_cg*p;
	  T fy = projectFunValue(y);
	  T fx = A.funcValue(result, B);
	  DEBUG_LOG("fx-fy: " << fx-fy);

	  result_is_prop = true;
	  while( result_is_prop && residual_out > tolerance_factor
			 && fy<=fx && iteration < max_iterations ){

		DEBUG_LOG(name << " cg_step");
		debug_fun(fun_val = fy;cout<<name<<setprecision(12)<<" func = "<<fun_val<<endl;);
		result = y;
		g -= alpha_cg*AP;
		projector.PHI(g, phi);  assert_ge(g.dot(phi),0);
		precond.solve(g,z, phi);   assert_eq(z, z);
		assert_ge_ext(g.dot(z),0,"\ng.dot(phi) = "<<g.dot(phi)<<"\n|phi|="<<phi.norm()<<"\n|beta|="<<beta.norm());
		const T beta = (z.dot(AP)) / (p.dot(AP)); assert_eq(beta, beta);
		p = z-beta*p;
		residual_out = computeGradients(g, false);
		DEBUG_LOG(name << " residual = " << residual_out);
		result_is_prop = proportional(result);

		A.multiply(p, AP);
		const T pd = p.dot(AP); assert_ge(pd, 0); // pd = p^t*A*p > 0
		alpha_cg = (z.dot(g))/pd; assert_ge(alpha_cg, 0);
		y = result-alpha_cg*p;
	    fx = fy;
	    fy = projectFunValue(y);

		num_cg ++;
		iteration ++;
	  }

	  return residual_out;
	}

  protected:
	//problem
	const MAT& A;
	const Vec& B;

	//parameter
	size_t max_iterations;
	T tolerance_factor;
	T gamma, alpha_bar;

	// internal structures
	PRECONDITION &precond;
	PROJECTOIN &projector;
	
	//temporary
	Vec g,p,z,beta,phi,gp,y;
	size_t iterations_out;
	T residual_out;

	// solving info
	int result_code;
	int num_cg, num_exp, numprop, iteration;
	T fun_val;
	std::string name;
  };

  // traditional MPRGP method framework.
  template <typename T, typename MAT, typename PROJECTOIN,typename PRECONDITION>
  class MPRGP:public MPRGPBase<T,MAT, PROJECTOIN, PRECONDITION>{

	typedef Eigen::Matrix<T,-1,1> Vec;
	typedef MPRGPBase<T,MAT,PROJECTOIN, PRECONDITION> MB;

  public:
	MPRGP(const MAT &A,const Vec &B,PRECONDITION &precond, 
		  PROJECTOIN &projector,const int max_it = 1000,const T tol=1e-3, Vec *ev = NULL):
	  MB(A,B,precond, projector, max_it, tol, ev){}

	int solve(Vec &result){

	  this->initialize(result);

	  for( MB::iteration = 0; MB::iteration < MB::max_iterations; MB::iteration++ ){

	  	DEBUG_LOG(MB::name << " step "<<MB::iteration);
	  	MB::residual_out = this->computeGradients(MB::g);

	  	if( MB::residual_out <= MB::tolerance_factor ){
	  	  MB::iterations_out = MB::iteration;
	  	  MB::result_code = 0;
	  	  break;
	  	}

	  	if( MB::proportional(result) ){

	  	  Vec& AP=MB::gp;
	  	  MB::A.multiply(MB::p,AP);
	  	  const T pd = MB::p.dot(AP); assert_ge(pd, 0); // pd = p^t*A*p > 0
	  	  const T alpha_cg = (MB::z.dot(MB::g)) / pd; assert_ge(alpha_cg, 0);
	  	  const T alpha_f = MB::projector.stepLimit(result,MB::p,alpha_cg); 
		  assert_ge(alpha_f, 0);

	  	  if(alpha_cg <= alpha_f){
 	  		this->CGStep(AP, alpha_cg, result);
	  	  }else{
	  		this->ExpStep(AP, alpha_f, result);
	  	  }
	  	}else{
	  	  this->PropStep(result);
	  	}
	  }

  	  MB::result_code = MB::residual_out <= MB::tolerance_factor ? 0 : -1;
  	  ERROR_LOG_COND(MB::name << " is not convergent with "<< MB::iteration << 
					 " iterations."<<endl,MB::result_code >= 0);

	  this->printSolveInfo();
	  return MB::result_code;
	}

  };

  // traditional monotonic MPRGP method framework.
  template <typename T, typename MAT, typename PROJECTOIN,typename PRECONDITION>
  class MPRGPMonotonic:public MPRGPBase<T,MAT, PROJECTOIN, PRECONDITION>{

  	typedef Eigen::Matrix<T,-1,1> Vec;
  	typedef MPRGPBase<T,MAT,PROJECTOIN, PRECONDITION> MB;

  public:
  	MPRGPMonotonic(const MAT &A,const Vec &B,PRECONDITION &precond, 
  				   PROJECTOIN &projector,const int max_it = 1000,const T tol=1e-3, Vec *ev = NULL):
  	  MB(A,B,precond, projector, max_it, tol, ev){}

  	int solve(Vec &result){

  	  this->initialize(result);

  	  for( MB::iteration = 0; MB::iteration < MB::max_iterations; ){

  	  	MB::residual_out = this->computeGradients(MB::g);
		DEBUG_LOG(MB::name << " residual = " << MB::residual_out);
  	  	if( MB::residual_out <= MB::tolerance_factor )
  	  	  break;

		bool result_is_prop = MB::proportional(result);
  	  	if(result_is_prop){
  		  MB::residual_out = this->CGMonotonicStep(result, result_is_prop);
		}

  	  	if( MB::residual_out <= MB::tolerance_factor )
  	  	  break;

  		if ( result_is_prop || MB::y.size() <= 0 || !(MB::projector.isFeasible(MB::y)) ){
  		  this->ExpMonotonicStep(result);
		  MB::iteration ++;
  	  	}else if( !result_is_prop ){
  		  MB::A.multiply(result, MB::g);
  		  MB::g -= MB::B;
  	  	  this->PropStep(result);
		  MB::iteration ++;
  	  	}
  	  }

	  MB::z = result;
	  MB::projector.project(MB::z, result);
	  
  	  MB::result_code = MB::residual_out <= MB::tolerance_factor ? 0 : -1;
  	  MB::iterations_out = MB::iteration;

  	  ERROR_LOG_COND(MB::name << " is not convergent with "<< MB::iteration << 
					 " iterations."<<endl,MB::result_code >= 0);

	  debug_fun(MB::projector.DECIDE_FACE(result));
  	  this->printSolveInfo();
  	  return MB::result_code;
  	}

  };

  // a matrix providing matrix-vector production, rows, and diagonal elements.
  template <typename T>
  class FixedSparseMatrix{

  public:
	FixedSparseMatrix(const SparseMatrix<T> &M):A(M){
	  assert_eq(A.rows(),A.cols());
	  getDiagonal(A, diag_A);
	}
	
	// result = A*x
	template <typename VEC,typename VEC_OUT>
	void multiply(const VEC& x,VEC_OUT& result)const{
	  assert_eq(A.cols(), x.size());
	  result = A*x;
	}
	
	// return 1/2*x^t*A*x-x^t*b
	template <typename VEC, typename VEC_SECOND>
	T funcValue(const VEC& x,const VEC_SECOND&b)const{
	  assert_eq(A.rows(),x.size());
	  assert_eq(A.cols(),x.size());
	  assert_eq(b.size(),x.size());
	  return 0.5f*x.dot(A*x)-x.dot(b);
	}

	int rows()const{
	  return A.rows();
	}

	double diag(const int i)const{
	  assert_in(i,0,diag_A.size()-1);
	  return diag_A[i];
	}

	const SparseMatrix<T> &getMatrix()const{
	  return A;
	}

  protected:
	SparseMatrix<T> A;
	Matrix<T,-1,1> diag_A;
  };

  // MPRGP solver for lower bound constraints: x >= L.
  template<typename T=double>
  class MPRGPLowerBound{

	typedef Eigen::Matrix<T,-1,1> Vec;
	
  public:
	template <typename MAT>
	static int solve(const MAT&A,const Vec&B,const Vec&L,Vec &x,const T tol=1e-3,const int max_it=1000){
	  
	  LowerBoundProjector<T> projector(L);
	  DiagonalInFacePreconSolver<T,MAT> precond(A, projector.getFace());
	  MPRGPMonotonic<T, MAT, LowerBoundProjector<T>, DiagonalInFacePreconSolver<T,MAT> > solver(A, B, precond, projector, max_it, tol);
	  return solver.solve(x);
	}
  };

  // MPRGP solver for box constraints: H >= x >= L.
  template<typename T=double>
  class MPRGPBoxBound{

	typedef Eigen::Matrix<T,-1,1> Vec;
	
  public:
	template <typename MAT> 
	static int solve(const MAT &A,const Vec &B, const Vec &L, const Vec &U, Vec &x, const T tol=1e-3, const int max_it = 1000){

	  BoxBoundProjector<T> projector(L, U);
	  DiagonalInFacePreconSolver<T,MAT> precond(A, projector.getFace());
	  MPRGP<T, MAT, BoxBoundProjector<T>, DiagonalInFacePreconSolver<T,MAT> > solver(A, B, precond, projector, max_it, tol);
	  const int rlst_code = solver.solve(x);
	  return rlst_code;
	}
  };

  // MPRGP solver for decoupled collision constraints: Jx >= c and J^t J should be diagonal.
  template<typename T=double>
  class MPRGPDecoupledCon{

	typedef Eigen::Matrix<T,-1,1> Vec;
	
  public:
	template <typename MAT, bool preconditioned = true>
	static int solve(const MAT &A,const Vec &B, DecoupledConProjector<T> &projector, Vec &x, 
					 const T tol=1e-3, const int max_it = 1000, 
					 const std::string solver_name = "MPRGP", Vec *ev = NULL){

	  assert_eq(A.rows(),B.size());
	  assert_eq(A.rows(),x.size());
	  typedef DiagonalDecouplePreconSolver<T,MAT, !preconditioned> Preconditioner;
	  Preconditioner precond(A, projector.getFace(), projector.getConMatrix());

	  typedef MPRGPMonotonic<T, MAT, DecoupledConProjector<T>, Preconditioner > MPRGPSolver;
	  MPRGPSolver solver(A, B, precond, projector, max_it, tol, ev);
	  solver.setName(solver_name);
	  const int rlst_code = solver.solve(x);
	  return rlst_code;
	}

	template <typename MAT, bool preconditioned = true>
	static int solve(const MAT &A,const Vec &B,const SparseMatrix<T> &J,const Vec &c,Vec &x, 
					 const T tol=1e-3, const int max_it = 1000, 
					 const std::string solver_name = "MPRGP", Vec *ev = NULL){

	  const SparseMatrix<T> JJt_mat = J*J.transpose();
	  assert_eq_ext(JJt_mat.nonZeros(), J.rows(), "Matrix J is not decoupled.\n" << J);
	  Vec JJt;
	  getDiagonal(JJt_mat, JJt);
	  DecoupledConProjector<T> projector(J, JJt, c);
	  const Vec init_x = x;
	  projector.project(init_x, x);
	  return solve<MAT, preconditioned>(A, B, projector, x, tol, max_it, solver_name, ev);
	}

	template <bool preconditioned = true>
	static int solve(const std::string file_name, Vec&x, 
					 const T tol=1e-3, const int max_it = 1000,
					 const std::string solver_name = "MPRGP", Vec *ev = NULL){

	  SparseMatrix<T> A, J;
	  Vec B, c;
	  int code = -1;
	  if ( loadQP(A, B, J, c, x, file_name) ){
		FixedSparseMatrix<T> FA(A);
		code = solve<FixedSparseMatrix<T>,preconditioned>(FA, B, J, c, x, tol, max_it, solver_name, ev);
	  }
	  return code;
	}
  };
  
}//end of namespace

#endif /* _MPRGPSOLVER_H_ */
