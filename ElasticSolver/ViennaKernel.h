#ifndef VIENNA_KERNEL_H
#define VIENNA_KERNEL_H

#include "solvers/MatVec.h"
#define VIENNACL_WITH_OPENCL 1
#define VIENNACL_WITH_EIGEN 1
#include <viennacl/matrix.hpp>
#include <viennacl/vector.hpp>
#include <viennacl/linalg/prod.hpp>
#include <viennacl/linalg/maxmin.hpp>
#include <viennacl/linalg/norm_2.hpp>
#include <viennacl/linalg/norm_1.hpp>
#include <viennacl/linalg/norm_inf.hpp>
#include <viennacl/linalg/inner_prod.hpp>
#include <viennacl/linalg/vector_operations.hpp>
#include <boost/algorithm/string.hpp>

PRJ_BEGIN

template <typename T>
struct ViennaKernel {
    template<typename T2>
    struct Rebind {
        typedef ViennaKernel<T2> value;
    };
    typedef T Scalar;
    typedef viennacl::vector<T> Vec;
    typedef viennacl::matrix<T,viennacl::column_major> Mat;
    static FORCE_INLINE T dot(const Vec& x, const Vec& y,sizeType n=-1) {
        return viennacl::linalg::inner_prod(x,y);
    }
    static FORCE_INLINE T norm(const Vec& x,sizeType n=-1) {
        return viennacl::linalg::norm_2(x);
    }
    static FORCE_INLINE T norm1(const Vec& x,sizeType n=-1) {
        return viennacl::linalg::norm_1(x);
    }
    static FORCE_INLINE void weightedPow(const Vec& x,const Vec& w,const T& p,Vec& y,sizeType n=-1) {
        if(!viennacl::ocl::current_context().has_program("weightedPowProgram")) {
            std::string weightedPowProgram=
                "#pragma OPENCL EXTENSION cl_khr_fp64: enable;\n"
                "__kernel void weightedPow(\n"
                "__global const SCALAR_TYPE *x,\n"
                "__global const SCALAR_TYPE *w,\n"
                "__global SCALAR_TYPE *y,\n"
                "unsigned int size,SCALAR_TYPE p) \n"
                "{\n"
                "	for(unsigned int i=get_global_id(0);\n"
                "		i<size;i+=get_global_size(0))\n"
                "	y[i]=w[i]*pow(fmax(x[i],(SCALAR_TYPE)1E-3f),p);\n"
                "};\n";
            boost::replace_all(weightedPowProgram,"SCALAR_TYPE",sizeof(T) == sizeof(double) ? "double" : "float");
            viennacl::ocl::current_context().add_program(weightedPowProgram,"weightedPowProgram");
        }
        viennacl::ocl::kernel& ker=
            viennacl::ocl::current_context().get_kernel("weightedPowProgram","weightedPow");
        viennacl::ocl::enqueue(ker(x,w,y,cl_uint(x.size()),p));
    }
    static FORCE_INLINE void copy(const Vec& x,Vec& y,sizeType n=-1) {
        if(y.size() != x.size())
            y.resize(x.size());
        y=x;
    }
    static FORCE_INLINE void copy(const typename Kernel<T>::Vec& x,Vec& y,sizeType n=-1) {
        if(y.size() != x.size())
            y.resize(x.size());
        viennacl::copy(x,y);
    }
    static FORCE_INLINE void copy(const Vec& x,typename Kernel<T>::Vec& y,sizeType n=-1) {
        if(y.size() != x.size())
            y.resize(x.size());
        viennacl::copy(x.begin(),x.end(),y.data());
    }
    static FORCE_INLINE void copy(const typename Kernel<T>::Mat& x,Mat& y) {
        if(y.size1() != x.rows() || y.size2() != x.cols())
            y.resize(x.rows(),x.cols());
        viennacl::copy(x,y);
    }
    static FORCE_INLINE void ncopy(const Vec& x,Vec& y,sizeType n=-1) {
        if(y.size() != x.size())
            y.resize(x.size());
        y=-x;
    }
    static FORCE_INLINE sizeType indexAbsMax(const Vec& x,sizeType n=-1) {
        Kernel<T>::Vec tmp(x.size());
        viennacl::copy(x.begin(),x.end(),tmp.data());
        return Kernel<T>::indexAbsMax(tmp);
    }
    static FORCE_INLINE T absMax(const Vec& x,sizeType n=-1) {
        return viennacl::linalg::norm_inf(x);
    }
    static FORCE_INLINE void scale(T alpha,Vec& y,sizeType n=-1) {
        y*=alpha;
    }
    static FORCE_INLINE void addScaled(const T& alpha, const Vec& x, Vec& y,sizeType n=-1) {
        y+=x*alpha;
    }
    static FORCE_INLINE void add(const Vec& x,const Vec& y,Vec& result,sizeType n=-1) {
        result=x+y;
    }
    static FORCE_INLINE void sub(const Vec& x,const Vec& y,Vec& result,sizeType n=-1) {
        result=x-y;
    }
    static FORCE_INLINE void cwiseProd(const Vec& x, const Vec& y, Vec& z,sizeType n=-1) {
        z=viennacl::linalg::element_prod(x,y);
    }
    static FORCE_INLINE void cwiseClamp(const T& l, const T& h, Vec& x,sizeType n=-1) {
        if(!viennacl::ocl::current_context().has_program("cwiseClampProgram")) {
            std::string cwiseClampProgram=
                "#pragma OPENCL EXTENSION cl_khr_fp64: enable;\n"
                "__kernel void cwiseClamp(\n"
                "__global SCALAR_TYPE *v,\n"
                "unsigned int size,\n"
                "SCALAR_TYPE l,SCALAR_TYPE h) \n"
                "{\n"
                "	for(unsigned int i=get_global_id(0);\n"
                "		i<size;i+=get_global_size(0))\n"
                "	v[i]=fmin(fmax(v[i],l),h);\n"
                "};\n";
            boost::replace_all(cwiseClampProgram,"SCALAR_TYPE",sizeof(T) == sizeof(double) ? "double" : "float");
            viennacl::ocl::current_context().add_program(cwiseClampProgram,"cwiseClampProgram");
        }
        viennacl::ocl::kernel& ker=
            viennacl::ocl::current_context().get_kernel("cwiseClampProgram","cwiseClamp");
        viennacl::ocl::enqueue(ker(x,cl_uint(x.size()),l,h));
    }
    static FORCE_INLINE void zero(Vec& v,sizeType n=-1) {
        v.clear();
    }
    static FORCE_INLINE void set(Vec& v,const T& val,sizeType n=-1) {
        std::vector<T> tmp(v.size(),val);
        viennacl::copy(tmp.begin(),tmp.end(),v.begin());
    }
    static FORCE_INLINE void mul(const Mat& A,const Vec& x,Vec& out) {
        out=viennacl::linalg::prod(A,x);
    }
    static FORCE_INLINE void mulT(const Mat& A,const Vec& x,Vec& out) {
        out=viennacl::linalg::prod(viennacl::trans(A),x);
    }
    static FORCE_INLINE sizeType rows(const Mat& A) {
        return A.size1();
    }
    static FORCE_INLINE sizeType cols(const Mat& A) {
        return A.size2();
    }
    static FORCE_INLINE T colNorm(const Mat& A,sizeType i) {
        Vec x;
        x=viennacl::column(A,(unsigned int)i);
        return norm(x);
    }
    template <typename A,typename B>
    static FORCE_INLINE void addMulT(A& x,const T& b,const B& c) {
        x+=b*c;
    }
    template <typename A,typename B>
    static FORCE_INLINE void subMulT(A& x,const T& b,const B& c) {
        x-=b*c;
    }
    static FORCE_INLINE T absV(const T& val) {
        return std::abs(val);
    }
    static FORCE_INLINE T ZeroV() {
        return (T)0.0f;
    }
};

PRJ_END

#endif