#ifndef FEM_UTIL_H
#define FEM_UTIL_H

#include "MathBasic.h"
#include "IO.h"
#include <Eigen/Sparse>
#include <boost/unordered_set.hpp>

PRJ_BEGIN

//hash
struct Hash {
    size_t operator()(const Vec2i& key) const;
    size_t operator()(const Vec3i& key) const;
    size_t operator()(const Vec4i& key) const;
};
//deformation gradient
static FORCE_INLINE void calcGComp2D(Eigen::Matrix<scalarD,9,12>& GComp,const Mat3d& d)
{
    //F=DH*invDM
    GComp.setZero();
    for(sizeType r=0; r<2; r++)
        for(sizeType c=0; c<2; c++)
            for(sizeType i=0; i<2; i++) {
                //DH(r,i)*invDM(i,c)=(x(i+1)(r)-x1(r))*invDM(i,c)
                GComp(c*3+r,(i+1)*3+r)+=d(i,c);
                GComp(c*3+r,r)-=d(i,c);
            }
}
static FORCE_INLINE void calcGComp3D(Eigen::Matrix<scalarD,9,12>& GComp,const Mat3d& d)
{
    //F=DH*invDM
    GComp.setZero();
    for(sizeType r=0; r<3; r++)
        for(sizeType c=0; c<3; c++)
            for(sizeType i=0; i<3; i++) {
                //DH(r,i)*invDM(i,c)=(x(i+1)(r)-x1(r))*invDM(i,c)
                GComp(c*3+r,(i+1)*3+r)+=d(i,c);
                GComp(c*3+r,r)-=d(i,c);
            }
}
//matrix
template <typename T>
static void addI3x3(vector<Eigen::Triplet<scalarD,sizeType> >& H,sizeType r,sizeType c,T coef)
{
    H.push_back(Eigen::Triplet<scalarD,sizeType>(r+0,c+0,coef));
    H.push_back(Eigen::Triplet<scalarD,sizeType>(r+1,c+1,coef));
    H.push_back(Eigen::Triplet<scalarD,sizeType>(r+2,c+2,coef));
}
template <typename MT>
static void addI(vector<Eigen::Triplet<scalarD,sizeType> >& H,sizeType r,sizeType c,const MT& coef)
{
    sizeType nr=coef.size();
    for(sizeType i=0; i<nr; i++)
        H.push_back(Eigen::Triplet<scalarD,sizeType>(r+i,c+i,coef[i]));
}
template <typename MT>
static void addBlock(vector<Eigen::Triplet<scalarD,sizeType> >& H,sizeType r,sizeType c,const MT& coef)
{
    sizeType nrR=coef.rows();
    sizeType nrC=coef.cols();
    for(sizeType i=0; i<nrR; i++)
        for(sizeType j=0; j<nrC; j++)
            H.push_back(Eigen::Triplet<scalarD,sizeType>(r+i,c+j,coef(i,j)));
}
static void addSparseBlock(vector<Eigen::Triplet<scalarD,sizeType> >& H,sizeType r,sizeType c,const Eigen::SparseMatrix<scalarD,0,sizeType>& mat)
{
    for(sizeType k=0; k<mat.outerSize(); ++k)
        for(Eigen::SparseMatrix<scalarD,0,sizeType>::InnerIterator it(mat,k); it; ++it)
            H.push_back(Eigen::Triplet<scalarD,sizeType>(it.row()+r,it.col()+c,it.value()));
}
template <typename T>
static void extend(T& m,sizeType r,sizeType c)
{
    T tmp=m;
    m.resize(m.rows()+r,m.cols()+c);
    m.block(0,0,tmp.rows(),tmp.cols())=tmp;
}
//mesh segmentation
struct FEMBody;
struct SegmentBody {
    static void segment(const FEMBody& body,sizeType& nrP,vector<sizeType>& keyMap,bool debugPt=false);
    static void writeVTK(const std::string& path,const FEMBody& body,const vector<sizeType>& keyMap);
    static void writeSelVTK(const std::string& path,const FEMBody& body,const vector<sizeType>& keyMap);
};
//nonnegative least square
struct NNLS {
    static bool solve(const Matd& A,const Cold& b,Cold& x);
};
//basis related
struct EigenSolver {
    static void solveEigen(const FEMBody& body,const boost::unordered_set<sizeType>& fixSet,Eigen::Matrix<scalarD,-1,1>& lambda,bool excludeRigid);
    static void makeOrthogonal(const Eigen::SparseMatrix<scalarD,0,sizeType>* A,Matd& U);
    static void scaleDiagonal(Eigen::SparseMatrix<scalarD,0,sizeType>& A,scalarD coef,scalarD minV);
};

PRJ_END

#endif
