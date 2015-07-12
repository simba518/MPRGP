#ifndef _CUFIXEDMAT_H_
#define _CUFIXEDMAT_H_
/**
 * @file   CuMatrix.h
 * @author simba <simba@simba>
 * @date   Sat Jul  2 16:43:56 2011
 * 
 * @brief  this is an simple matrix class that allocate memory from device memory.\n
 * used for cuda blas functions.
 * @see CuVector.h
*/
#include <assert.h>
#include <cuda_runtime_api.h>
#include <stdio.h>

#include <iostream>
using namespace std;

#include <CuFixedVec.h>

namespace LSW_CUDA_MATH{
  
  template <typename real>
  class CuFixedMat{
	
  public:
	//*_data should be in device memory
	CuFixedMat(int _r,int _c,real *_data){
	  this->r = _r;
	  this->c = _c;
	  this->data = _data;
	}

	bool sameDim(const CuFixedMat<real> &other)const{
	  return (this->rows() == other.rows() && this->columns() == other.columns());
	}

	int rows()const{ return r;}

	int cols()const{ return c;}

	real *Data(){ return data;}

	const real *Data()const{ return data;}

	const real &element(int i,int j)const{
	  assert (i>= 0 && j>= 0&&i < rows() && j < columns ());
	  int index = i + j*rows();
	  return data[index];
	}

	real &element(int i,int j){
	  assert (i>= 0 && j>= 0&&i < rows() && j < columns ());
	  int index = i + j*rows();
	  return data[index];
	}

	/// @brief initialize from host memory datas.
	bool copyFromHost(const real *h_data,int _r,int _c){

	  if(this->rows() == _r && this->columns() == _c){
		return ( cudaMemcpy(data,h_data,_r*_c*sizeof(real),cudaMemcpyHostToDevice) == cudaSuccess );
	  }
	  return false;
	}

	/// @brief copy datas to  host memory.
	/// @note h_data should be allocated outside.
	bool copyToHost(real *h_data)const{
	  if( h_data != NULL && this->Data() != NULL){
		return ( cudaMemcpy(h_data,data,r*c*sizeof(real),cudaMemcpyDeviceToHost) == cudaSuccess );
	  }
	  return false;
	}

	~CuFixedMat(){
	  //data should not be set to NULL
	  //because it should be released in the child class
	}
	
	CuFixedMat<real> & operator = (CuFixedMat<real> &other){
	  
	  this->r = other.rows();
	  this->c = other.columns();
	  this->data = other.Data();
	  return *this;
	}
	///return column[col_id] as an vector,
	///@note the return CuVector should not be resize or release.
	CuFixedVec<real> operator [] (int col_id){
	  assert (col_id >=0 && col_id < c);
	  return CuFixedVec<real>(r,&(this->element(0,col_id)));
	}
	const CuFixedVec<real> operator [] (int col_id)const{
	  assert (col_id >=0 && col_id < c);
	  return CuFixedVec<real>(r,&(this->element(0,col_id)));
	}
	
	///return the sub matrix cosisted by columns between [begin,begin+size-1].
	const CuFixedMat<real> columns(int begin,int size)const{
	  
	  assert(begin>=0 && size>0 && begin+size<=c);
	  return CuFixedMat<real> (rows(),size,&(data[begin*rows()]));
	}
	CuFixedMat<real> columns(int begin,int size){
	  
	  assert(begin>=0 && size>0 && begin+size<=c);
	  return CuFixedMat<real> (rows(),size,&(data[begin*rows()]));
	}

	///*this += other*s
	CuFixedMat<real> & add(const CuFixedMat<real> &other,real s=1.0f){
	  assert(this->sameDim(other));
	  int data_len = rows()*columns();
	  if(data_len > 0){
		bool succ = axpy(other.Data(),this->Data(),data_len,s);
		assert(succ);
	  }
	  return *this;
	}
	CuFixedMat<real> & operator += (const CuFixedMat<real> &other){
	  return add(other);
	}
	CuFixedMat<real> & operator -= (const CuFixedMat<real> &other){
	  return add(other,(real)(-1.0f));
	}
	///*this *= s
	CuFixedMat<real> & operator *= (const real s){
	  int size = rows()*columns();
	  if(size > 0){
		bool succ = scal(this->Data(),size,s);
		assert(succ);
	  }
	  return *this;
	}

	bool setAll0s(){
	  bool succ = true;
	  int data_len = rows()*cols();
	  if(data_len>0){
		succ = ( cudaMemset(this->data,0,data_len*sizeof(real)) == cudaSuccess );
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
	  printf("r = %d, c = %d",r,c);
	  if(r*c > 0){
		real *tempt_m = new real[r*c];
		this->copyToHost(tempt_m);
		for (int i=0; i< r; i++){
		  printf("\n");
		  for (int j=0; j < c; ++j){
			printf(format,tempt_m[i+j*r]);
		  }
		}
		printf("\n");
		delete [] tempt_m;
	  }
	}

	friend ostream & operator<<(ostream& theStream,const CuFixedMat<real> &me){
	  me.print();
	  return theStream;
	}


  protected:
	int r;
	int c;
	real *data;
	
  };

  typedef CuFixedMat<double> CuFixedMatXd;
  typedef CuFixedMat<float> CuFixedMatXf;

}//end of namespace

#endif /*_CUFIXEDMAT_H_*/
