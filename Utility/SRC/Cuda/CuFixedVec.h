#ifndef _CUFIXEDVECTOR_H_
#define _CUFIXEDVECTOR_H_

/**
 * @file   CuFixedVec.h
 * @author simba <simba@simba>
 * @date   Sat Jul  2 17:07:06 2011
 * 
 * @brief  it is an fixed size vector, fixed size means 
 * the length of the vector can't be changed. it is used 
 * for visiting the segment vector or column vector of 
 * CuVector or CuMatrix.
 * @see CuMatrix.h and CuVector.h
 */

#include <assert.h>
#include <cuda_runtime_api.h>
#include <stdio.h>

#include <iostream>
using namespace std;

#include "WarpedCuBlas.h"

namespace LSW_CUDA_MATH{
  
  template <typename real>
  class CuFixedVec{
		
  public:
	CuFixedVec(int len,real *_data){
	  length = len;
	  data = _data;
	}

	int size()const{ return length;}

	real *Data(){ return data;}

	const real *Data()const{ return data;}

	const real &operator [] (int i)const{
	  assert (i>= 0 && i < size());
	  return data[i];
	}

	real &operator [] (int i){
	  assert (i>= 0 && i < size());
	  return data[i];
	}

	/// @brief initialize from host memory datas.
	bool copyFromHost(const real *h_data,int len){

	  if(len <= length){
		return ( cudaMemcpy(data,h_data,len*sizeof(real),cudaMemcpyHostToDevice) == cudaSuccess );
	  }
	  return false;
	}

	/// @brief copy datas to  host memory.
	/// @note h_data should be allocated outside.
	bool copyToHost(real *h_data)const{
	  if( h_data != NULL && this->Data() != NULL){
		return ( cudaMemcpy(h_data,data,length*sizeof(real),cudaMemcpyDeviceToHost) == cudaSuccess );
	  }
	  return false;
	}

	///the datas in "other" will copy to this vector
	CuFixedVec<real> &copy(const CuFixedVec<real> &other){
	  
	  assert(this->size() == other.size());
	  bool succ = ( cudaMemcpy(this->data,other.Data(),length*sizeof(real),cudaMemcpyDeviceToDevice) == cudaSuccess );
	  assert (succ);
	  return (*this);
	}

	///the datas won't be copied, they shall the same data.
	CuFixedVec<real> & operator = (CuFixedVec<real> &other){
	  this->length = other.size();
	  this->data = other.Data();
	}

	const CuFixedVec<real> segment(int start,int len)const{
	  assert (start>=0 && len >0 && start+len <= this->size());
	  return CuFixedVec<real>(len,&(data[start]));
	}

	CuFixedVec<real> segment(int start,int len){
	  assert (start>=0 && len >0 && start+len <= this->size());
	  return CuFixedVec<real>(len,&(data[start]));
	}

	CuFixedVec<real> head(int len)const{
	  assert (len < this->size() && len > 0);
	  return segment(0,len);
	}

	CuFixedVec<real> head(int len){
	  assert (len < this->size() && len > 0);
	  return segment(0,len);	  
	}

	CuFixedVec<real> tail(int len)const{
	  assert (len < this->size() && len > 0);
	  return segment(this->size()-len,len);	  
	}

	CuFixedVec<real> tail(int len){
	  assert (len < this->size() && len > 0);
	  return segment(this->size()-len,len);	  	  
	}

	///*this += other*s
	CuFixedVec<real> & add(const CuFixedVec<real> &other,real s=1.0f){
	  assert(other.size() == this->size());
	  if(other.size() > 0){
		bool succ = axpy(other.Data(),this->Data(),other.size(),s);
		assert(succ);
	  }
	  return *this;
	}
	CuFixedVec<real> & operator += (const CuFixedVec<real> &other){
	  return add(other);
	}
	CuFixedVec<real> & operator -= (const CuFixedVec<real> &other){
	  return add(other,(real)(-1.0f));
	}
	///*this *= s
	CuFixedVec<real> & operator *= (const real s){
	  return this->scale(s);
	}
	CuFixedVec<real> & scale(const real s){

	  if(this->size() > 0){
		bool succ = scal(this->Data(),this->size(),s);
		assert(succ);
	  }
	  return *this;
	}
	
	

	///(*this).dot(v)
	///@note the return  value is on the host memory.
	real dot_host(const CuFixedVec<real> &v)const{
	  assert(this->size() == v.size());
	  real h_reslt = 0;
	  if(v.size()>0){
		real *d_reslt = NULL;
		bool succ = (cudaMalloc((void**)&d_reslt,sizeof(real)) == cudaSuccess);
		assert(succ);
		succ = gemm(data,1,size(),v.Data(),v.size(),1,d_reslt,'n','n',(real)1.0f,(real)0.0f);
		assert(succ);
		succ = ( cudaMemcpy(&h_reslt,d_reslt,sizeof(real),cudaMemcpyDeviceToHost) == cudaSuccess );
		assert(succ);
	  }
	  return h_reslt;
	}

	real norm_host()const{
	  return dot_host(*this);
	}

	bool setAll0s(){
	  bool succ = true;
	  if(this->size()>0){
		succ = ( cudaMemset(this->data,0,length*sizeof(real)) == cudaSuccess );
	  }
	  return succ;
	}

	bool setZero(){
	  return setAll0s();
	}

	void print(const char *name=NULL,const char *format="%4.4f\t")const{
	  
	  if(name != NULL){
		printf("%s\n",name);
	  }
	  printf("len = %d\n",length);
	  if(length > 0){
		real *tempt_v = new real[length];
		this->copyToHost(tempt_v);
		for (int i=0; i< length; i++){
		  printf(format,tempt_v[i]);
		}
		printf("\n");
		delete [] tempt_v;
	  }
	}

	friend ostream & operator<<(ostream& theStream,const CuFixedVec<real> &me){

	  me.print();
	  return theStream;
	}

	~CuFixedVec(){
	  //data should not be set to NULL
	  //because it should be released in the child class
	}

  protected:
	int length;
	real *data;
  };

  typedef CuFixedVec<double> CuFixedVecXd;
  typedef CuFixedVec<float> CuFixedVecXf;
  
}//end of namespace

#endif /*_CUFIXEDVECTOR_H_*/
