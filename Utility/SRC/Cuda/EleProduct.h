#ifndef _ELEPRODUCT_H_
#define _ELEPRODUCT_H_

#include <stdio.h>

namespace LSW_CUDA_MATH{

  class EleProduct{
	
  public:
	static EleProduct *getInstance(){
	  
	  if(p_instance == NULL){
		p_instance = new EleProduct();
	  }
	  return p_instance;
	}
	~EleProduct(){
	  
	  if(p_instance != NULL){
		delete p_instance;
	  }
	  p_instance = NULL;
	} 

	bool product(const double *d_v1,const double *d_v2,double *d_v1v2,int len)const;
	bool product(const float *d_v1,const float *d_v2,float *d_v1v2,int len)const;

  protected:
	EleProduct(){}
	
  private:
	static EleProduct *p_instance;
  };
    
}//end of namespace

#endif /*_ELEPRODUCT_H_*/
