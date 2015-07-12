#include <stdio.h>
#include <cuda_runtime_api.h>

namespace LSW_CUDA{

  template<typename real>
  __global__ void EleProductKernelFun(const real* v1, const real* v2, real* v1v2,int N){

	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < N){
	  v1v2[i] = v1[i]*v2[i];
	}
  }

  /** 
   * use cuda to compute the per element product of vector v1 and v2, 
   * store the result in v1v2.
   * v1v2[i] = v1[i]*v2[i].
   * It is used to compute the matrix-vector dot product, when the 
   * matrix is an diaginal matrix.
   *
   * @param d_v1 vector on device
   * @param d_v2 vector on device
   * @param d_v1v2 result of the product, on device, allocated outside.
   * @param len the length of the vectors.(all have the same length).
   * 
   * @return true if parameters are valid and compute success.
   */
  template<typename real>
  bool CuEleProduct(const real *d_v1,const real *d_v2,real *d_v1v2,int len){

	if(len == 0){
	  return true;
	}

	//check parameters
	bool succ = (d_v1!=NULL && d_v1v2!=NULL && d_v2!=NULL && len > 0);
	if (!succ){
	  printf("\nerror: CuEleProduct(...) invalid parameters!\n");
	  return false;
	}

	int N = len;
	int threadsPerBlock = 256;
	int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;

	//running kernel function on the device
	EleProductKernelFun<<<blocksPerGrid, threadsPerBlock>>>(d_v1, d_v2, d_v1v2, N);

	//check errors
	cudaError_t err = cudaGetLastError();
	if( cudaSuccess != err) {
	  succ = false;
	  printf("\nerror: CuEleProduct(...) kernel function EleProductKernelFun(..) failed!\n");
	  printf("cuda error:%s",cudaGetErrorString(err));
	}
	return succ;
  }

  extern"C"
  bool EleProductD(const double *d_v1,const double *d_v2,double *d_v1v2,int len){
  	return CuEleProduct(d_v1,d_v2,d_v1v2,len);
  }

  extern"C"
  bool EleProductF(const float *d_v1,const float *d_v2,float *d_v1v2,int len){
  	return CuEleProduct(d_v1,d_v2,d_v1v2,len);
  }

}//end of namespace
