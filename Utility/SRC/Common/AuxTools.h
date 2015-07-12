#ifndef _AUXTOOLS_H_
#define _AUXTOOLS_H_

#include <fstream>
#include <string>
#include <boost/lexical_cast.hpp>
#include <Log.h>

namespace UTILITY{

  /////////////////////convertion///////////////////////////
#define TOSTR(input) boost::lexical_cast<std::string>(input)
#define TOINT(input) boost::lexical_cast<int>(input)
#define TOFLOAT(input) boost::lexical_cast<double>(input)
  ////////////////////loop,size/////////////////////////////
#define SZ(X) ((int)(X).size())
#define ALL(X) (X).begin(), (X).end()

#define REP(I, N) for (int I = 0; I < (N); ++I)
#define REPP(I, A, B) for (int I = (A); I < (B); ++I)
#define REPC(I, C) for (int I = 0; !(C); ++I)

  ////////////////////io///////////////////////////////////
#define OUTFILE(file,name)												\
	std::ofstream file;													\
	file.open(std::string(name).c_str());								\
	ERROR_LOG_COND("failed to open file for output: "<<name,file.is_open()); 

#define INFILE(file,name)												\
	std::ifstream file;													\
	file.open(std::string(name).c_str());								\
	ERROR_LOG_COND("failed to open file for readin: "<<name,file.is_open()); 

  template <class SCALAR>
  inline void printVec(const SCALAR *vec,const int len,const std::string space=","){
	REP (i,len){
	  std::cout << vec[i];
	  (i != len-1) ? std::cout<<space:std::cout << std::endl;
	}
  }

  template <class VECTOR>
  inline void printVec(const VECTOR &vec,int len,const std::string space=","){
	if(len > 0) printVec(&vec[0],len,space);
  }

  template <class VECTOR>
  inline void printVec(const VECTOR &vec,const std::string space=","){
	printVec(vec,vec.size(),space);
  }

  ////////////////////math/////////////////////////////////
  template <class SCALAR>
  inline SCALAR norm2(const SCALAR *vec,const size_t len){
	SCALAR sum = 0.0f;
	REP(i,len) sum += (vec[i]*vec[i]);
	return sqrt(sum);
  }

  template <class VECTOR>
  inline double norm2(const VECTOR &vec,const size_t len){
	return (len>0 ? norm2(&vec[0],len):0);
  }

  template <class VECTOR>
  inline double norm2(const VECTOR &vec){
	return norm2(vec,vec.size());
  }
}
#endif /* _AUXTOOLS_H_ */
