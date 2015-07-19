#ifndef _SIMPLEQP_H_
#define _SIMPLEQP_H_

#include <Eigen/Dense>
#include <assertext.h>
#include "QuadProg.h"
using namespace Eigen;

typedef int64_t sizeType;
typedef Eigen::Matrix<double,-1,-1> Matd;
typedef Eigen::Matrix<double,-1, 1> Cold;

namespace MATH{

  //min  0.5 x^T A x - b^T x, s.t. Jx >= c0
  template<class MATRIX_A>
  void solveQP(const MATRIX_A& A,const Cold& b,const Matd& J,const Cold& c0,Cold& x){

	assert_eq(A.rows(), A.cols());
	assert_eq(b.size(), A.rows());
	assert_eq(J.cols(), b.size());
	assert_eq(c0.size(), J.rows());

	sizeType N=A.rows();
	sizeType C=J.rows();
	QuadProgPP::Matrix<double> G(N,N);
	QuadProgPP::Vector<double> g0(N);
	QuadProgPP::Matrix<double> CE(N,0);
	QuadProgPP::Vector<double> ce0(0);
	QuadProgPP::Matrix<double> CI(N,C);
	QuadProgPP::Vector<double> ci0(C);
	for(sizeType r=0;r<N;r++){
	  g0[r]=-b[r];
	  for(sizeType c=0;c<N;c++)
		G[r][c]=A(r,c);
	  for(sizeType c=0;c<C;c++)
		{
		  CI[r][c]=J(c,r);
		  ci0[c]=-c0[c];
		}
	}

	QuadProgPP::Vector<double> x0(N);
	QuadProgPP::QuadProg::solve_quadprog(G,g0,CE,ce0,CI,ci0,x0);
	x.resize(A.rows());
	for(sizeType r=0;r<N;r++)
	  x[r]=x0[r];
  }

  template<class MATRIX_A>
  void checkKKT(const Matd& A,const Cold& b,const Matd& J,const Cold& c0,const Cold& x){

	//solve dual
	sizeType C=J.rows();

	Matd invA=A.inverse();
	Matd AL=J*invA*J.transpose();
	Cold bL=c0-J*invA*b;
	Matd JL=Matd::Identity(C,C);
	Cold c0L=Cold::Zero(C);

	Cold lambda=Cold::Zero(C);
	solveQP(AL,bL,JL,c0L,lambda);
	Cold xTest=invA*(J.transpose()*lambda+b);

	//test
	std::cout << "Lambda" << std::endl;
	std::cout << lambda << std::endl << std::endl;
	
	std::cout << "Constraint Error" << std::endl;
	std::cout << (J*xTest-c0) << std::endl << std::endl;

	std::cout << "X Error" << std::endl;
	std::cout << (x-xTest).norm() << std::endl << std::endl;
  }

  // int main(){

  // 	sizeType nrN=100;
  // 	sizeType nrC=10;

  // 	Matd tmp=Matd::Random(nrN,nrN);
  // 	Matd A=tmp.transpose()*tmp;
  // 	Cold b=Cold::Random(nrN);

  // 	Matd J=Matd::Random(nrC,nrN);
  // 	Cold c0=Cold::Random(nrC);

  // 	Cold x=Cold::Zero(nrN);
  // 	solveQP(A,b,J,c0,x);
  // 	checkKKT(A,b,J,c0,x);
  // }

}//end of namespace

#endif /* _SIMPLEQP_H_ */
