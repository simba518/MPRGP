#ifndef KRYLOV_MATRIX_H
#define KRYLOV_MATRIX_H

#include "MatVec.h"

PRJ_BEGIN

template <typename T,typename KERNEL_TYPE=Kernel<T> >
struct KrylovMatrix {
public:
    typedef typename KERNEL_TYPE::Vec Vec;
    virtual ~KrylovMatrix() {}
    virtual Vec operator*(const Vec& b) const {
        Vec tmp(b.size());
        multiply(b,tmp);
        return tmp;
    }
    virtual void multiply(const Vec& b,Vec& out) const=0;
    virtual sizeType n() const =0;
};
template <typename T,typename KERNEL_TYPE=Kernel<T> >
struct IdentityKrylovMatrix : public KrylovMatrix<T,KERNEL_TYPE> {
public:
    typedef typename KERNEL_TYPE::Vec Vec;
	IdentityKrylovMatrix(sizeType n):_n(n){}
    virtual void multiply(const Vec& b,Vec& out) const {out=b;}
    virtual sizeType n() const {return _n;}
	sizeType _n;
};
template <typename T,typename KERNEL_TYPE=Kernel<T> >
struct DefaultKrylovMatrix : public KrylovMatrix<T,KERNEL_TYPE> {
public:
    typedef typename KERNEL_TYPE::Vec Vec;
    DefaultKrylovMatrix(const FixedSparseMatrix<T,KERNEL_TYPE> &matrix) {
        _fixedMatrix=matrix;
    }
    DefaultKrylovMatrix(const SparseMatrix<T,KERNEL_TYPE> &matrix) {
        _fixedMatrix.constructFromMatrix(matrix);
    }
	DefaultKrylovMatrix(const sizeType nr,std::vector<Eigen::Triplet<T,sizeType> >& trips){
		_fixedMatrix.resize(nr);
		_fixedMatrix.buildFromTriplets(trips.begin(),trips.end());
	}
    virtual void multiply(const Vec& b,Vec& out) const {
        _fixedMatrix.multiply(b,out);
    }
    virtual sizeType n() const {
        return _fixedMatrix.rows();
    }
    FixedSparseMatrix<T,KERNEL_TYPE> _fixedMatrix;
};
template <typename T,typename KERNEL_TYPE=Kernel<T> >
struct EigenKrylovMatrix : public KrylovMatrix<T,KERNEL_TYPE> {
public:
    typedef typename KERNEL_TYPE::Vec Vec;
	EigenKrylovMatrix(){}
    EigenKrylovMatrix(const Eigen::SparseMatrix<T,0,sizeType> &matrix):_fixedMatrix(matrix){}
    virtual void multiply(const Vec& b,Vec& out) const {
        out=_fixedMatrix*b;
    }
    virtual sizeType n() const {
        return _fixedMatrix.rows();
    }
	Eigen::SparseMatrix<T,0,sizeType> _fixedMatrix;
};
template <typename T,typename KERNEL_TYPE=Kernel<T> >
struct DenseEigenKrylovMatrix : public KrylovMatrix<T,KERNEL_TYPE> {
public:
    typedef typename KERNEL_TYPE::Vec Vec;
	DenseEigenKrylovMatrix(){}
    DenseEigenKrylovMatrix(const Eigen::Matrix<T,-1,-1> &matrix):_fixedMatrix(matrix){}
    virtual void multiply(const Vec& b,Vec& out) const {
        out=_fixedMatrix*b;
    }
    virtual sizeType n() const {
        return _fixedMatrix.rows();
    }
    Eigen::Matrix<T,-1,-1> _fixedMatrix;
};
template <typename T,typename KERNEL_TYPE=Kernel<T> >
struct DefaultEigenKrylovMatrix : public KrylovMatrix<T,KERNEL_TYPE> {
public:
    typedef typename KERNEL_TYPE::Vec Vec;
    DefaultEigenKrylovMatrix(const Eigen::SparseMatrix<T,0,sizeType> &matrix):_fixedMatrix(matrix){}
    virtual void multiply(const Vec& b,Vec& out) const {
        out=_fixedMatrix*b;
    }
    virtual sizeType n() const {
        return _fixedMatrix.rows();
    }
    const Eigen::SparseMatrix<T,0,sizeType>& _fixedMatrix;
};
template <typename T,typename KERNEL_TYPE=Kernel<T> >
struct DefaultDenseEigenKrylovMatrix : public KrylovMatrix<T,KERNEL_TYPE> {
public:
    typedef typename KERNEL_TYPE::Vec Vec;
    DefaultDenseEigenKrylovMatrix(const Eigen::Matrix<T,-1,-1> &matrix):_fixedMatrix(matrix){}
    virtual void multiply(const Vec& b,Vec& out) const {
        out=_fixedMatrix*b;
    }
    virtual sizeType n() const {
        return _fixedMatrix.rows();
    }
    const Eigen::Matrix<T,-1,-1>& _fixedMatrix;
};
template <typename T,typename KERNEL_TYPE=Kernel<T> >
struct DiagonalEigenKrylovMatrix : public KrylovMatrix<T,KERNEL_TYPE> {
public:
    typedef typename KERNEL_TYPE::Vec Vec;
	DiagonalEigenKrylovMatrix(){}
    DiagonalEigenKrylovMatrix(const Eigen::SparseMatrix<T,0,sizeType> &matrix):_fixedMatrix(matrix){}
    virtual void multiply(const Vec& b,Vec& out) const {
        out=_fixedMatrix*b;
    }
    virtual sizeType n() const {
        return _fixedMatrix.rows();
    }
	Eigen::DiagonalMatrix<T,-1> _fixedMatrix;
};

PRJ_END

#endif