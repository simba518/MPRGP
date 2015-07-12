#include "EleProduct.h"

using namespace LSW_CUDA_MATH;

EleProduct *EleProduct::p_instance = NULL;

extern"C"
bool EleProductD(const double *d_v1,const double *d_v2,double *d_v1v2,int len);

extern"C"
bool EleProductF(const float *d_v1,const float *d_v2,float *d_v1v2,int len);

bool EleProduct::product(const double *d_v1,const double *d_v2,double *d_v1v2,int len)const{
  
  return EleProductD(d_v1,d_v2,d_v1v2,len);
}

bool EleProduct::product(const float *d_v1,const float *d_v2,float *d_v1v2,int len)const{

  return EleProductF(d_v1,d_v2,d_v1v2,len);
}
