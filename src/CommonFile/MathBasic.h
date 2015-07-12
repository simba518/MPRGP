#ifndef MATH_BASIC_H
#define MATH_BASIC_H

#include "Config.h"
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <omp.h>
#include <float.h>
#include <Eigen/Eigen>
#include <stdio.h>
#include <sys/time.h>
#include <boost/timer/timer.hpp>
#include <boost/pool/pool_alloc.hpp>

#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

PRJ_BEGIN

using namespace std;

//types
#ifdef M_PI
#undef M_PI
#endif

template <typename T>
struct ScalarUtil;
template <>
struct ScalarUtil<float> {
    typedef Eigen::Vector4f ScalarVec4;
    typedef Eigen::Vector3f ScalarVec3;
    typedef Eigen::Vector2f ScalarVec2;

    typedef Eigen::Matrix2f ScalarMat2;
    typedef Eigen::Matrix3f ScalarMat3;
    typedef Eigen::Matrix4f ScalarMat4;

    static float scalar_max;
    static float scalar_eps;
};
template <>
struct ScalarUtil<double> {
    typedef Eigen::Vector4d ScalarVec4;
    typedef Eigen::Vector3d ScalarVec3;
    typedef Eigen::Vector2d ScalarVec2;

    typedef Eigen::Matrix2d ScalarMat2;
    typedef Eigen::Matrix3d ScalarMat3;
    typedef Eigen::Matrix4d ScalarMat4;

    static double scalar_max;
    static double scalar_eps;
};
template <>
struct ScalarUtil<sizeType> {
    typedef Eigen::Matrix<sizeType,4,1> ScalarVec4;
    typedef Eigen::Matrix<sizeType,3,1> ScalarVec3;
    typedef Eigen::Matrix<sizeType,2,1> ScalarVec2;

    typedef Eigen::Matrix<sizeType,2,2> ScalarMat2;
    typedef Eigen::Matrix<sizeType,3,3> ScalarMat3;
    typedef Eigen::Matrix<sizeType,4,4> ScalarMat4;

    static sizeType scalar_max;
    static sizeType scalar_eps;
};

template<typename T> sizeType getEigenStride()
{
    vector<T,Eigen::aligned_allocator<T> > tester(2);
    return (sizeType)(((char*)&(tester[1]))-((char*)&(tester[0])));
}

#ifndef DOUBLE_PRECISION
#define SCALAR_MAX FLT_MAX
#define SCALAR_EPS FLT_EPSILON
#define M_PI 3.14159265358979323846f
#define EPS 1E-5f
typedef float scalar;
typedef Eigen::Vector4f Vec4;
typedef Eigen::Vector3f Vec3;
typedef Eigen::Vector2f Vec2;
typedef Eigen::VectorXf Vecx;

typedef Eigen::Matrix2f Mat2;
typedef Eigen::Matrix3f Mat3;
typedef Eigen::Matrix4f Mat4;

typedef Eigen::Quaternionf Quat;
#else
#define SCALAR_MAX DBL_MAX
#define SCALAR_EPS DBL_EPSILON
#define M_PI 3.14159265358979323846
#define EPS 1E-9
typedef double scalar;
typedef Eigen::Vector4d Vec4;
typedef Eigen::Vector3d Vec3;
typedef Eigen::Vector2d Vec2;
typedef Eigen::VectorXd Vecx;

typedef Eigen::Matrix2d Mat2;
typedef Eigen::Matrix3d Mat3;
typedef Eigen::Matrix4d Mat4;

typedef Eigen::Quaterniond Quat;
#endif

typedef double scalarD;
typedef float scalarF;

typedef Eigen::Matrix2d Mat2d;
typedef Eigen::Matrix3d Mat3d;
typedef Eigen::Matrix4d Mat4d;
typedef Eigen::Quaterniond Quatd;

typedef Eigen::Matrix<scalarD,-1, 1> Cold;
typedef Eigen::Matrix<scalarD, 1,-1> Rowd;
typedef Eigen::Matrix<scalarD,-1,-1> Matd;

typedef Eigen::Matrix2f Mat2f;
typedef Eigen::Matrix3f Mat3f;
typedef Eigen::Matrix4f Mat4f;
typedef Eigen::Quaternionf Quatf;

typedef Eigen::Matrix<scalarF,-1, 1> Colf;
typedef Eigen::Matrix<scalarF, 1,-1> Rowf;
typedef Eigen::Matrix<scalarF,-1,-1> Matf;

typedef Eigen::Vector4d Vec4d;
typedef Eigen::Vector3d Vec3d;
typedef Eigen::Vector2d Vec2d;

typedef Eigen::Vector4f Vec4f;
typedef Eigen::Vector3f Vec3f;
typedef Eigen::Vector2f Vec2f;

typedef Eigen::Matrix<sizeType, 4, 1> Vec4i;
typedef Eigen::Matrix<sizeType, 3, 1> Vec3i;
typedef Eigen::Matrix<sizeType, 2, 1> Vec2i;

typedef Eigen::Matrix<sizeType,-1, 1> Coli;
typedef Eigen::Matrix<sizeType, 1,-1> Rowi;
typedef Eigen::Matrix<sizeType,-1,-1> Mati;

typedef Eigen::Matrix<unsigned char, 4, 1> Vec4uc;
typedef Eigen::Matrix<unsigned char, 3, 1> Vec3uc;
typedef Eigen::Matrix<unsigned char, 2, 1> Vec2uc;

typedef Eigen::Matrix<char, 4, 1> Vec4c;
typedef Eigen::Matrix<char, 3, 1> Vec3c;
typedef Eigen::Matrix<char, 2, 1> Vec2c;

FORCE_INLINE bool compL(const Vec3i& a,const Vec3i& b)
{
    return a.x() < b.x() && a.y() < b.y() && a.z() < b.z();
}
FORCE_INLINE bool compLE(const Vec3i& a,const Vec3i& b)
{
    return a.x() <= b.x() && a.y() <= b.y() && a.z() <= b.z();
}
FORCE_INLINE bool compG(const Vec3i& a,const Vec3i& b)
{
    return a.x() > b.x() && a.y() > b.y() && a.z() > b.z();
}
FORCE_INLINE bool compGE(const Vec3i& a,const Vec3i& b)
{
    return a.x() >= b.x() && a.y() >= b.y() && a.z() >= b.z();
}
FORCE_INLINE bool compL(const Vec2i& a,const Vec2i& b)
{
    return a.x() < b.x() && a.y() < b.y();
}
FORCE_INLINE bool compLE(const Vec2i& a,const Vec2i& b)
{
    return a.x() <= b.x() && a.y() <= b.y();
}
FORCE_INLINE bool compG(const Vec2i& a,const Vec2i& b)
{
    return a.x() > b.x() && a.y() > b.y();
}
FORCE_INLINE bool compGE(const Vec2i& a,const Vec2i& b)
{
    return a.x() >= b.x() && a.y() >= b.y();
}

FORCE_INLINE bool compL(const Vec3f& a,const Vec3f& b)
{
    return a.x() < b.x() && a.y() < b.y() && a.z() < b.z();
}
FORCE_INLINE bool compLE(const Vec3f& a,const Vec3f& b)
{
    return a.x() <= b.x() && a.y() <= b.y() && a.z() <= b.z();
}
FORCE_INLINE bool compG(const Vec3f& a,const Vec3f& b)
{
    return a.x() > b.x() && a.y() > b.y() && a.z() > b.z();
}
FORCE_INLINE bool compGE(const Vec3f& a,const Vec3f& b)
{
    return a.x() >= b.x() && a.y() >= b.y() && a.z() >= b.z();
}
FORCE_INLINE bool compL(const Vec2f& a,const Vec2f& b)
{
    return a.x() < b.x() && a.y() < b.y();
}
FORCE_INLINE bool compLE(const Vec2f& a,const Vec2f& b)
{
    return a.x() <= b.x() && a.y() <= b.y();
}
FORCE_INLINE bool compG(const Vec2f& a,const Vec2f& b)
{
    return a.x() > b.x() && a.y() > b.y();
}
FORCE_INLINE bool compGE(const Vec2f& a,const Vec2f& b)
{
    return a.x() >= b.x() && a.y() >= b.y();
}

FORCE_INLINE bool compL(const Vec3d& a,const Vec3d& b)
{
    return a.x() < b.x() && a.y() < b.y() && a.z() < b.z();
}
FORCE_INLINE bool compLE(const Vec3d& a,const Vec3d& b)
{
    return a.x() <= b.x() && a.y() <= b.y() && a.z() <= b.z();
}
FORCE_INLINE bool compG(const Vec3d& a,const Vec3d& b)
{
    return a.x() > b.x() && a.y() > b.y() && a.z() > b.z();
}
FORCE_INLINE bool compGE(const Vec3d& a,const Vec3d& b)
{
    return a.x() >= b.x() && a.y() >= b.y() && a.z() >= b.z();
}
FORCE_INLINE bool compL(const Vec2d& a,const Vec2d& b)
{
    return a.x() < b.x() && a.y() < b.y();
}
FORCE_INLINE bool compLE(const Vec2d& a,const Vec2d& b)
{
    return a.x() <= b.x() && a.y() <= b.y();
}
FORCE_INLINE bool compG(const Vec2d& a,const Vec2d& b)
{
    return a.x() > b.x() && a.y() > b.y();
}
FORCE_INLINE bool compGE(const Vec2d& a,const Vec2d& b)
{
    return a.x() >= b.x() && a.y() >= b.y();
}

FORCE_INLINE scalarF compMin(const scalarF& a,const scalarF& b){return min(a,b);}
FORCE_INLINE scalarF compMax(const scalarF& a,const scalarF& b){return max(a,b);}
FORCE_INLINE scalarD compMin(const scalarD& a,const scalarD& b){return min(a,b);}
FORCE_INLINE scalarD compMax(const scalarD& a,const scalarD& b){return max(a,b);}

FORCE_INLINE Vec2i compMin(const Vec2i& a,const Vec2i& b)
{
    return Vec2i(min(a.x(),b.x()),min(a.y(),b.y()));
}
FORCE_INLINE Vec3i compMin(const Vec3i& a,const Vec3i& b)
{
    return Vec3i(min(a.x(),b.x()),min(a.y(),b.y()),min(a.z(),b.z()));
}
FORCE_INLINE Vec4i compMin(const Vec4i& a,const Vec4i& b)
{
    return Vec4i(min(a.x(),b.x()),min(a.y(),b.y()),min(a.z(),b.z()),min(a.w(),b.w()));
}

FORCE_INLINE Vec2i compMax(const Vec2i& a,const Vec2i& b)
{
    return Vec2i(max(a.x(),b.x()),max(a.y(),b.y()));
}
FORCE_INLINE Vec3i compMax(const Vec3i& a,const Vec3i& b)
{
    return Vec3i(max(a.x(),b.x()),max(a.y(),b.y()),max(a.z(),b.z()));
}
FORCE_INLINE Vec4i compMax(const Vec4i& a,const Vec4i& b)
{
    return Vec4i(max(a.x(),b.x()),max(a.y(),b.y()),max(a.z(),b.z()),max(a.w(),b.w()));
}

FORCE_INLINE Vec2f compMin(const Vec2f& a,const Vec2f& b)
{
    return Vec2f(min(a.x(),b.x()),min(a.y(),b.y()));
}
FORCE_INLINE Vec3f compMin(const Vec3f& a,const Vec3f& b)
{
    return Vec3f(min(a.x(),b.x()),min(a.y(),b.y()),min(a.z(),b.z()));
}
FORCE_INLINE Vec4f compMin(const Vec4f& a,const Vec4f& b)
{
    return Vec4f(min(a.x(),b.x()),min(a.y(),b.y()),min(a.z(),b.z()),min(a.w(),b.w()));
}

FORCE_INLINE Vec2f compMax(const Vec2f& a,const Vec2f& b)
{
    return Vec2f(max(a.x(),b.x()),max(a.y(),b.y()));
}
FORCE_INLINE Vec3f compMax(const Vec3f& a,const Vec3f& b)
{
    return Vec3f(max(a.x(),b.x()),max(a.y(),b.y()),max(a.z(),b.z()));
}
FORCE_INLINE Vec4f compMax(const Vec4f& a,const Vec4f& b)
{
    return Vec4f(max(a.x(),b.x()),max(a.y(),b.y()),max(a.z(),b.z()),max(a.w(),b.w()));
}

FORCE_INLINE Vec2uc compMin(const Vec2uc& a,const Vec2uc& b)
{
    return Vec2uc(min(a.x(),b.x()),min(a.y(),b.y()));
}
FORCE_INLINE Vec3uc compMin(const Vec3uc& a,const Vec3uc& b)
{
    return Vec3uc(min(a.x(),b.x()),min(a.y(),b.y()),min(a.z(),b.z()));
}
FORCE_INLINE Vec4uc compMin(const Vec4uc& a,const Vec4uc& b)
{
    return Vec4uc(min(a.x(),b.x()),min(a.y(),b.y()),min(a.z(),b.z()),min(a.w(),b.w()));
}

FORCE_INLINE Vec2uc compMax(const Vec2uc& a,const Vec2uc& b)
{
    return Vec2uc(max(a.x(),b.x()),max(a.y(),b.y()));
}
FORCE_INLINE Vec3uc compMax(const Vec3uc& a,const Vec3uc& b)
{
    return Vec3uc(max(a.x(),b.x()),max(a.y(),b.y()),max(a.z(),b.z()));
}
FORCE_INLINE Vec4uc compMax(const Vec4uc& a,const Vec4uc& b)
{
    return Vec4uc(max(a.x(),b.x()),max(a.y(),b.y()),max(a.z(),b.z()),max(a.w(),b.w()));
}

FORCE_INLINE Vec2c compMin(const Vec2c& a,const Vec2c& b)
{
    return Vec2c(min(a.x(),b.x()),min(a.y(),b.y()));
}
FORCE_INLINE Vec3c compMin(const Vec3c& a,const Vec3c& b)
{
    return Vec3c(min(a.x(),b.x()),min(a.y(),b.y()),min(a.z(),b.z()));
}
FORCE_INLINE Vec4c compMin(const Vec4c& a,const Vec4c& b)
{
    return Vec4c(min(a.x(),b.x()),min(a.y(),b.y()),min(a.z(),b.z()),min(a.w(),b.w()));
}

FORCE_INLINE Vec2c compMax(const Vec2c& a,const Vec2c& b)
{
    return Vec2c(max(a.x(),b.x()),max(a.y(),b.y()));
}
FORCE_INLINE Vec3c compMax(const Vec3c& a,const Vec3c& b)
{
    return Vec3c(max(a.x(),b.x()),max(a.y(),b.y()),max(a.z(),b.z()));
}
FORCE_INLINE Vec4c compMax(const Vec4c& a,const Vec4c& b)
{
    return Vec4c(max(a.x(),b.x()),max(a.y(),b.y()),max(a.z(),b.z()),max(a.w(),b.w()));
}

FORCE_INLINE Vec2i ceil(const Vec2f& vec)
{
    return Vec2i((sizeType)std::ceil(vec.x()),(sizeType)std::ceil(vec.y()));
}
FORCE_INLINE Vec3i ceil(const Vec3f& vec)
{
    return Vec3i((sizeType)std::ceil(vec.x()),(sizeType)std::ceil(vec.y()),(sizeType)std::ceil(vec.z()));
}
FORCE_INLINE Vec4i ceil(const Vec4f& vec)
{
    return Vec4i((sizeType)std::ceil(vec.x()),(sizeType)std::ceil(vec.y()),(sizeType)std::ceil(vec.z()),(sizeType)std::ceil(vec.w()));
}

FORCE_INLINE Vec2i floor(const Vec2f& vec)
{
    return Vec2i((sizeType)std::floor(vec.x()),(sizeType)std::floor(vec.y()));
}
FORCE_INLINE Vec3i floor(const Vec3f& vec)
{
    return Vec3i((sizeType)std::floor(vec.x()),(sizeType)std::floor(vec.y()),(sizeType)std::floor(vec.z()));
}
FORCE_INLINE Vec4i floor(const Vec4f& vec)
{
    return Vec4i((sizeType)std::floor(vec.x()),(sizeType)std::floor(vec.y()),(sizeType)std::floor(vec.z()),(sizeType)std::floor(vec.w()));
}

FORCE_INLINE Vec2i round(const Vec2f& vec)
{
    return Vec2i((sizeType)std::floor(vec.x()+0.5f),(sizeType)std::floor(vec.y()+0.5f));
}
FORCE_INLINE Vec3i round(const Vec3f& vec)
{
    return Vec3i((sizeType)std::floor(vec.x()+0.5f),(sizeType)std::floor(vec.y()+0.5f),(sizeType)std::floor(vec.z()+0.5f));
}
FORCE_INLINE Vec4i round(const Vec4f& vec)
{
    return Vec4i((sizeType)std::floor(vec.x()+0.5f),(sizeType)std::floor(vec.y()+0.5f),(sizeType)std::floor(vec.z()+0.5f),(sizeType)std::floor(vec.w()+0.5f));
}

FORCE_INLINE Vec2d compMin(const Vec2d& a,const Vec2d& b)
{
    return Vec2d(min(a.x(),b.x()),min(a.y(),b.y()));
}
FORCE_INLINE Vec3d compMin(const Vec3d& a,const Vec3d& b)
{
    return Vec3d(min(a.x(),b.x()),min(a.y(),b.y()),min(a.z(),b.z()));
}
FORCE_INLINE Vec4d compMin(const Vec4d& a,const Vec4d& b)
{
    return Vec4d(min(a.x(),b.x()),min(a.y(),b.y()),min(a.z(),b.z()),min(a.w(),b.w()));
}

FORCE_INLINE Vec2d compMax(const Vec2d& a,const Vec2d& b)
{
    return Vec2d(max(a.x(),b.x()),max(a.y(),b.y()));
}
FORCE_INLINE Vec3d compMax(const Vec3d& a,const Vec3d& b)
{
    return Vec3d(max(a.x(),b.x()),max(a.y(),b.y()),max(a.z(),b.z()));
}
FORCE_INLINE Vec4d compMax(const Vec4d& a,const Vec4d& b)
{
    return Vec4d(max(a.x(),b.x()),max(a.y(),b.y()),max(a.z(),b.z()),max(a.w(),b.w()));
}

FORCE_INLINE Vec2i ceil(const Vec2d& vec)
{
    return Vec2i((sizeType)std::ceil(vec.x()),(sizeType)std::ceil(vec.y()));
}
FORCE_INLINE Vec3i ceil(const Vec3d& vec)
{
    return Vec3i((sizeType)std::ceil(vec.x()),(sizeType)std::ceil(vec.y()),(sizeType)std::ceil(vec.z()));
}
FORCE_INLINE Vec4i ceil(const Vec4d& vec)
{
    return Vec4i((sizeType)std::ceil(vec.x()),(sizeType)std::ceil(vec.y()),(sizeType)std::ceil(vec.z()),(sizeType)std::ceil(vec.w()));
}

FORCE_INLINE Vec2i floor(const Vec2d& vec)
{
    return Vec2i((sizeType)std::floor(vec.x()),(sizeType)std::floor(vec.y()));
}
FORCE_INLINE Vec3i floor(const Vec3d& vec)
{
    return Vec3i((sizeType)std::floor(vec.x()),(sizeType)std::floor(vec.y()),(sizeType)std::floor(vec.z()));
}
FORCE_INLINE Vec4i floor(const Vec4d& vec)
{
    return Vec4i((sizeType)std::floor(vec.x()),(sizeType)std::floor(vec.y()),(sizeType)std::floor(vec.z()),(sizeType)std::floor(vec.w()));
}

FORCE_INLINE Vec2i round(const Vec2d& vec)
{
    return Vec2i((sizeType)std::floor(vec.x()+0.5),(sizeType)std::floor(vec.y()+0.5));
}
FORCE_INLINE Vec3i round(const Vec3d& vec)
{
    return Vec3i((sizeType)std::floor(vec.x()+0.5),(sizeType)std::floor(vec.y()+0.5),(sizeType)std::floor(vec.z()+0.5));
}
FORCE_INLINE Vec4i round(const Vec4d& vec)
{
    return Vec4i((sizeType)std::floor(vec.x()+0.5),(sizeType)std::floor(vec.y()+0.5),(sizeType)std::floor(vec.z()),(sizeType)std::floor(vec.w()+0.5));
}

template <typename T,int DIM=3>
struct BBox {
    typedef Eigen::Matrix<T,DIM,1> PT;
    typedef Eigen::Matrix<T,2,1> PT2;
    BBox() {
        reset();
    }
    virtual ~BBox() {}
    BBox(const PT& p):_minC(p),_maxC(p) {}
    BBox(const PT& minC,const PT& maxC):_minC(minC),_maxC(maxC) {}
    template <typename T2>
    BBox(const BBox<T2,DIM>& other) {
        copy(other);
    }
    static BBox createMM(const PT& minC,const PT& maxC) {
        return BBox(minC,maxC);
    }
    static BBox createME(const PT& minC,const PT& extent) {
        return BBox(minC,minC+extent);
    }
    static BBox createCE(const PT& center,const PT& extent) {
        return BBox(center-extent,center+extent);
    }
    BBox getIntersect(const BBox& other) const {
        return createMM(compMax(_minC,other._minC),compMin(_maxC,other._maxC));
    }
    BBox getUnion(const BBox& other) const {
        return createMM(compMin(_minC,other._minC),compMax(_maxC,other._maxC));
    }
    BBox getUnion(const PT& point) const {
        return createMM(compMin(_minC,point),compMax(_maxC,point));
    }
    BBox getUnion(const PT& ctr,const T& rad) const {
        return createMM(compMin(_minC,ctr-PT::Constant(rad)),compMax(_maxC,ctr+PT::Constant(rad)));
    }
    void setIntersect(const BBox& other) {
        _minC=compMax(_minC,other._minC);
        _maxC=compMin(_maxC,other._maxC);
    }
    void setUnion(const BBox& other) {
        _minC=compMin(_minC,other._minC);
        _maxC=compMax(_maxC,other._maxC);
    }
    void setUnion(const PT& point) {
        _minC=compMin(_minC,point);
        _maxC=compMax(_maxC,point);
    }
    void setUnion(const PT& ctr,const T& rad) {
        _minC=compMin(_minC,ctr-PT::Constant(rad));
        _maxC=compMax(_maxC,ctr+PT::Constant(rad));
    }
    void setPoints(const PT& a,const PT& b,const PT& c) {
        _minC=compMin(compMin(a,b),c);
        _maxC=compMax(compMax(a,b),c);
    }
	PT minCorner() const{return _minC;}
	PT maxCorner() const{return _maxC;}
    void enlargedEps(T eps) {
        PT d=(_maxC-_minC)*(eps*0.5f);
        _minC-=d;
        _maxC+=d;
    }
    BBox enlargeEps(T eps) const {
        PT d=(_maxC-_minC)*(eps*0.5f);
        return createMM(_minC-d,_maxC+d);
    }
    void enlarged(T len,const sizeType d=DIM) {
        for(sizeType i=0;i<d;i++){
            _minC[i]-=len;
            _maxC[i]+=len;
        }
    }
    BBox enlarge(T len,const sizeType d=DIM) const {
        BBox ret=createMM(_minC,_maxC);
        ret.enlarged(len,d);
        return ret;
    }
    PT lerp(const PT& frac) const {
        return _maxC*frac-_minC*(frac-1.0f);
    }
    bool empty() const {
        return !compL(_minC,_maxC);
    }
    template <int DIM2>
    bool containDim(const PT& point) const {
        for(int i=0;i<DIM2;i++)
            if(_minC[i] > point[i] || _maxC[i] < point[i])
                return false;
        return true;
    }
    bool contain(const BBox& other,const sizeType d=DIM) const {
        for(int i=0;i<d;i++)
            if(_minC[i] > other._minC[i] || _maxC[i] < other._maxC[i])
                return false;
        return true;
    }
    bool contain(const PT& point,const sizeType d=DIM) const {
        for(int i=0;i<d;i++)
            if(_minC[i] > point[i] || _maxC[i] < point[i])
                return false;
        return true;
    }
    bool contain(const PT& point,const T& rad,const sizeType d=DIM) const {
        for(int i=0;i<d;i++)
            if(_minC[i]+rad > point[i] || _maxC[i]-rad < point[i])
                return false;
        return true;
    }
    void reset() {
        _minC=PT::Constant( ScalarUtil<T>::scalar_max);
        _maxC=PT::Constant(-ScalarUtil<T>::scalar_max);
    }
    PT getExtent() const {
        return _maxC-_minC;
    }
	T distTo(const BBox& other,const sizeType d=DIM) const
	{
		PT dist=PT::Zero();
        for(sizeType i=0; i<d; i++) {
            if (other._maxC[i] < _minC[i])
                dist[i] = other._maxC[i] - _minC[i];
            else if (other._minC[i] > _maxC[i])
                dist[i] = other._minC[i] - _maxC[i];
        }
        return dist.norm();
	}
    T distTo(const PT& pt,const sizeType d=DIM) const {
        return sqrt(distToSqr(pt,d));
    }
    T distToSqr(const PT& pt,const sizeType d=DIM) const {
        PT dist=PT::Zero();
        for(sizeType i=0; i<d; i++) {
            if (pt[i] < _minC[i])
                dist[i] = pt[i] - _minC[i];
            else if (pt[i] > _maxC[i])
                dist[i] = pt[i] - _maxC[i];
        }
        return dist.squaredNorm();
    }
    PT closestTo(const PT& pt,const sizeType d=DIM) const {
        PT dist(pt);
        for(sizeType i=0; i<d; i++) {
            if (pt[i] < _minC[i])
                dist[i] = _minC[i];
            else if (pt[i] > _maxC[i])
                dist[i] = _maxC[i];
        }
        return dist;
    }
    bool intersect(const PT& p,const PT& q,const sizeType d=DIM) const {
        const T lo=1-5*ScalarUtil<T>::scalar_eps;
        const T hi=1+5*ScalarUtil<T>::scalar_eps;

        T s=0, t=1;
        for(sizeType i=0; i<d; ++i) {
            if(p[i]<q[i]) {
                T d=q[i]-p[i];
                T s0=(_minC[i]-p[i])/d, t0=(_maxC[i]-p[i])/d;
                if(s0>s) s=s0;
                if(t0<t) t=t0;
            } else if(p[i]>q[i]) {
                T d=q[i]-p[i];
                T s0=lo*(_maxC[i]-p[i])/d, t0=hi*(_minC[i]-p[i])/d;
                if(s0>s) s=s0;
                if(t0<t) t=t0;
            } else {
                if(p[i]<_minC[i] || p[i]>_maxC[i])
                    return false;
            }

            if(s>t)
                return false;
        }
        return true;
    }
    bool intersect(const PT& p,const PT& q,T& s,T& t,const sizeType d=DIM) const {
        const T lo=1-5*ScalarUtil<T>::scalar_eps;
        const T hi=1+5*ScalarUtil<T>::scalar_eps;

        s=0;
        t=1;
        for(sizeType i=0; i<d; ++i) {
            if(p[i]<q[i]) {
                T d=q[i]-p[i];
                T s0=lo*(_minC[i]-p[i])/d, t0=hi*(_maxC[i]-p[i])/d;
                if(s0>s) s=s0;
                if(t0<t) t=t0;
            } else if(p[i]>q[i]) {
                T d=q[i]-p[i];
                T s0=lo*(_maxC[i]-p[i])/d, t0=hi*(_minC[i]-p[i])/d;
                if(s0>s) s=s0;
                if(t0<t) t=t0;
            } else {
                if(p[i]<_minC[i] || p[i]>_maxC[i])
                    return false;
            }

            if(s>t)
                return false;
        }
        return true;
    }
    bool intersect(const BBox& other,const sizeType& d=DIM) const {
        for(sizeType i=0; i<d; i++) {
            if(_maxC[i] < other._minC[i] || other._maxC[i] < _minC[i])
                return false;
        }
        return true;
        //return compLE(_minC,other._maxC) && compLE(other._minC,_maxC);
    }
    PT2 project(const PT& a,const sizeType d=DIM) const {
        PT ctr=(_minC+_maxC)*0.5f;
        T ctrD=a.dot(ctr);
        T delta=0.0f;
        ctr=_maxC-ctr;
        for(sizeType i=0; i<d; i++)
            delta+=std::abs(ctr[i]*a[i]);
        return PT2(ctrD-delta,ctrD+delta);
    }
    template<typename T2>
    BBox& copy(const BBox<T2,DIM>& other) {
        for(sizeType i=0; i<DIM; i++) {
            _minC[i]=(T)other._minC[i];
            _maxC[i]=(T)other._maxC[i];
        }
        return *this;
    }
    T perimeter(const sizeType d=DIM) const
    {
        PT ext=getExtent();
        if(d <= 2)
            return ext.sum()*2.0f;
        else{ 
            ASSERT(d == 3);
            return (ext[0]*ext[1]+ext[1]*ext[2]+ext[0]*ext[2])*2.0f;
        }
    }
    PT _minC;
    PT _maxC;
};
typedef BBox<scalarF> BBoxf;
typedef BBox<scalarD> BBoxd;

template <typename T> FORCE_INLINE T mmin(const T& a,const T& b)
{
    return a < b ? a : b;
}
template <typename T> FORCE_INLINE T mmax(const T& a,const T& b)
{
    return a > b ? a : b;
}

struct ErfCalc {
    static const sizeType ncof=28;
    static const scalarD cof[28];
    scalarD erf(scalarD x) const {
        if (x >=0.0) return 1.0 - erfccheb(x);
        else return erfccheb(-x) - 1.0;
    }
    scalarD erfc(scalarD x) const {
        if (x >= 0.0) return erfccheb(x);
        else return 2.0 - erfccheb(-x);
    }
    scalarD erfccheb(scalarD z) const {
        sizeType j;
        scalarD t,ty,tmp,d=0.0,dd=0.0;
        if (z < 0.0)
            throw("erfccheb requires nonnegative argument");
        t = 2.0/(2.0+z);
        ty = 4.0*t - 2.0;
        for (j=ncof-1; j>0; j--) {
            tmp = d;
            d = ty*d - dd + cof[j];
            dd = tmp;
        }
        return t*exp(-z*z + 0.5*(cof[0] + ty*d) - dd);
    }
    scalarD inverfc(scalarD p) const {
        scalarD x,err,t,pp;
        if (p >= 2.0) return -100.0;
        if (p <= 0.0) return 100.0;
        pp = (p < 1.0)? p : 2.0 - p;
        t = sqrt(-2.0*log(pp/2.0));
        x = -0.70711*((2.30753+t*0.27061)/(1.+t*(0.99229+t*0.04481)) - t);
        for (sizeType j=0; j<2; j++) {
            err = erfc(x) - pp;
            x += err/(1.12837916709551257*exp(-(x*x))-x*err);
        }
        return (p < 1.0? x : -x);
    }
    scalarD inverf(scalarD p) const {
        return inverfc(1.-p);
    }
};
FORCE_INLINE scalarD erf(const scalarD& val)
{
    return ErfCalc().erf(val);   //only double accepted
}
template <typename T>
FORCE_INLINE typename ScalarUtil<T>::ScalarVec3 randBary()
{
    typename ScalarUtil<T>::ScalarVec3 ret;

    ret.x()=rand()/(T)RAND_MAX;
    ret.y()=rand()/(T)RAND_MAX;
    ret.y()*=(1.0f-ret.x());
    ret.z()=1.0f-ret.x()-ret.y();

    return ret;
}
class RandHash
{
public:
    // returns repeatable stateless pseudo-random number in [0,1]
    static FORCE_INLINE scalarD randhashd(sizeType seed) {
        return randhashInner((unsigned int)seed)/(scalarD)UINT_MAX;
    }
    static FORCE_INLINE scalarF randhashf(sizeType seed) {
        return randhashInner((unsigned int)seed)/(scalarF)UINT_MAX;
    }
    static FORCE_INLINE scalar randhash(sizeType seed) {
        return randhashInner((unsigned int)seed)/(scalar)UINT_MAX;
    }
    static FORCE_INLINE sizeType randhashI(sizeType seed) {
        sizeType ret=0;
        if(sizeof(ret) == 8) {
            ret=(sizeType)randhashInner((unsigned int)seed);
            ret<<=32;
            ret+=(sizeType)randhashInner((unsigned int)seed);
        } else {
            ret=(sizeType)randhashInner((unsigned int)seed);
        }
        return ret;
    }
    // returns repeatable stateless pseudo-random number in [a,b]
    static FORCE_INLINE scalarD randhashd(sizeType seed, scalarD a, scalarD b) {
        return (b-a)*randhashInner((unsigned int)seed)/(scalarD)UINT_MAX + a;
    }
    static FORCE_INLINE scalarF randhashf(sizeType seed, scalarF a, scalarF b) {
        return (b-a)*randhashInner((unsigned int)seed)/(scalarF)UINT_MAX + a;
    }
    static FORCE_INLINE scalar randhash(sizeType seed, scalar a, scalar b) {
        return (b-a)*randhashInner((unsigned int)seed)/(scalar)UINT_MAX + a;
    }
protected:
    static FORCE_INLINE unsigned int randhashInner(unsigned int seed) {
        unsigned int i=(seed^0xA3C59AC3u)*2654435769u;
        i^=(i>>16);
        i*=2654435769u;
        i^=(i>>16);
        i*=2654435769u;
        return i;
    }
    // the inverse of randhash
    static FORCE_INLINE unsigned int unhashInner(unsigned int h) {
        h*=340573321u;
        h^=(h>>16);
        h*=340573321u;
        h^=(h>>16);
        h*=340573321u;
        h^=0xA3C59AC3u;
        return h;
    }
};
FORCE_INLINE sizeType sgn(const scalarD& val)
{
    return val > 0.0f ? 1 : val < 0.0f ? -1 : 0;
}
FORCE_INLINE sizeType sgn(const scalarF& val)
{
    return val > 0.0f ? 1 : val < 0.0f ? -1 : 0;
}
FORCE_INLINE bool isInf( const scalarD& x )
{
    return ( ( x - x ) != 0.0f );
}
FORCE_INLINE bool isInf( const scalarF& x )
{
    return ( ( x - x ) != 0.0f );
}
FORCE_INLINE bool isNan( const scalarD& x )
{
    return ( ( ! ( x < 0.0f ) ) && ( ! ( x >= 0.0f ) ) );
}
FORCE_INLINE bool isNan( const scalarF& x )
{
    return ( ( ! ( x < 0.0f ) ) && ( ! ( x >= 0.0f ) ) );
}
FORCE_INLINE sizeType makePOT(sizeType val,sizeType decimate)
{
    ASSERT(decimate > 1)

    sizeType ret=1;
    while(ret < val)
        ret*=decimate;
    return ret;
}
FORCE_INLINE sizeType mod(sizeType a,sizeType b,const scalar &bInv)
{
    sizeType n=sizeType(a*bInv);
    a-=n*b;
    if(a<0)
        a+=b;
    return a;
}
FORCE_INLINE scalarD minmod( const scalarD& a,const scalarD& b,const scalarD& c)
{
    if(a > 0.0f && b > 0.0f && c > 0.0f)
        return min<scalarD>(min<scalarD>(a,b),c);
    else if(a < 0.0f && b < 0.0f && c < 0.0f)
        return max<scalarD>(max<scalarD>(a,b),c);
    else
        return 0.0f;
}
FORCE_INLINE scalarF minmod( const scalarF& a,const scalarF& b,const scalarF& c)
{
    if(a > 0.0f && b > 0.0f && c > 0.0f)
        return min<scalarF>(min<scalarF>(a,b),c);
    else if(a < 0.0f && b < 0.0f && c < 0.0f)
        return max<scalarF>(max<scalarF>(a,b),c);
    else
        return 0.0f;
}

template<typename RET>
RET countBits(const unsigned char& v)
{
    static const unsigned char BitsSetTable256[256] = {
#   define B2(n) n,     n+1,     n+1,     n+2
#   define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#   define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
        B6(0), B6(1), B6(1), B6(2)
    };
    return BitsSetTable256[v];
}
template<typename RET>
RET countBits(const sizeType& v)
{
    RET ret=0;
    unsigned char* vc=(unsigned char*)&v;
    for(sizeType i=0; i<sizeof(sizeType); i++)
        ret+=countBits<RET>(vc[i]);
    return ret;
}

FORCE_INLINE Vec2i makePOT(const Vec2i& val,sizeType decimate)
{
    return Vec2i(makePOT(val.x(),decimate),makePOT(val.y(),decimate));
}
FORCE_INLINE Vec3i makePOT(const Vec3i& val,sizeType decimate)
{
    return Vec3i(makePOT(val.x(),decimate),makePOT(val.y(),decimate),makePOT(val.z(),decimate));
}
FORCE_INLINE Vec4i makePOT(const Vec4i& val,sizeType decimate)
{
    return Vec4i(makePOT(val.x(),decimate),makePOT(val.y(),decimate),makePOT(val.z(),decimate),makePOT(val.w(),decimate));
}

FORCE_INLINE Vec2i operator>>(const Vec2i& in,const sizeType& nr)
{
    return Vec2i(in.x() >> nr,in.y() >> nr);
}
FORCE_INLINE Vec3i operator>>(const Vec3i& in,const sizeType& nr)
{
    return Vec3i(in.x() >> nr,in.y() >> nr,in.z() >> nr);
}
FORCE_INLINE Vec4i operator>>(const Vec4i& in,const sizeType& nr)
{
    return Vec4i(in.x() >> nr,in.y() >> nr,in.z() >> nr,in.w() >> nr);
}
FORCE_INLINE Vec2i operator&(const Vec2i& in,const sizeType& nr)
{
    return Vec2i(in.x()&nr,in.y()&nr);
}
FORCE_INLINE Vec3i operator&(const Vec3i& in,const sizeType& nr)
{
    return Vec3i(in.x()&nr,in.y()&nr,in.z()&nr);
}
FORCE_INLINE Vec4i operator&(const Vec4i& in,const sizeType& nr)
{
    return Vec4i(in.x()&nr,in.y()&nr,in.z()&nr,in.w()&nr);
}

FORCE_INLINE bool isInf( const Vec2f& x )
{
    return isInf(x.x()) || isInf(x.y());
}
FORCE_INLINE bool isNan( const Vec2f& x )
{
    return isNan(x.x()) || isNan(x.y());
}
FORCE_INLINE bool isInf( const Vec3f& x )
{
    return isInf(x.x()) || isInf(x.y()) || isInf(x.z());
}
FORCE_INLINE bool isNan( const Vec3f& x )
{
    return isNan(x.x()) || isNan(x.y()) || isNan(x.z());
}
FORCE_INLINE bool isInf( const Vec4f& x )
{
    return isInf(x.x()) || isInf(x.y()) || isInf(x.z()) || isInf(x.w());
}
FORCE_INLINE bool isNan( const Vec4f& x )
{
    return isNan(x.x()) || isNan(x.y()) || isNan(x.z()) || isNan(x.w());
}
FORCE_INLINE Vec2f minmodV(const Vec2f& a,const Vec2f& b,const Vec2f& c)
{
    return Vec2f(minmod(a.x(),b.x(),c.x()),
                 minmod(a.y(),b.y(),c.y()));
}
FORCE_INLINE Vec3f minmodV(const Vec3f& a,const Vec3f& b,const Vec3f& c)
{
    return Vec3f(minmod(a.x(),b.x(),c.x()),
                 minmod(a.y(),b.y(),c.y()),
                 minmod(a.z(),b.z(),c.z()));
}
FORCE_INLINE Vec4f minmodV(const Vec4f& a,const Vec4f& b,const Vec4f& c)
{
    return Vec4f(minmod(a.x(),b.x(),c.x()),
                 minmod(a.y(),b.y(),c.y()),
                 minmod(a.z(),b.z(),c.z()),
                 minmod(a.w(),b.w(),c.w()));
}

FORCE_INLINE bool isInf( const Vec2d& x )
{
    return isInf(x.x()) || isInf(x.y());
}
FORCE_INLINE bool isNan( const Vec2d& x )
{
    return isNan(x.x()) || isNan(x.y());
}
FORCE_INLINE bool isInf( const Vec3d& x )
{
    return isInf(x.x()) || isInf(x.y()) || isInf(x.z());
}
FORCE_INLINE bool isNan( const Vec3d& x )
{
    return isNan(x.x()) || isNan(x.y()) || isNan(x.z());
}
FORCE_INLINE bool isInf( const Vec4d& x )
{
    return isInf(x.x()) || isInf(x.y()) || isInf(x.z()) || isInf(x.w());
}
FORCE_INLINE bool isNan( const Vec4d& x )
{
    return isNan(x.x()) || isNan(x.y()) || isNan(x.z()) || isNan(x.w());
}
FORCE_INLINE Vec2d minmodV(const Vec2d& a,const Vec2d& b,const Vec2d& c)
{
    return Vec2d(minmod(a.x(),b.x(),c.x()),
                 minmod(a.y(),b.y(),c.y()));
}
FORCE_INLINE Vec3d minmodV(const Vec3d& a,const Vec3d& b,const Vec3d& c)
{
    return Vec3d(minmod(a.x(),b.x(),c.x()),
                 minmod(a.y(),b.y(),c.y()),
                 minmod(a.z(),b.z(),c.z()));
}
FORCE_INLINE Vec4d minmodV(const Vec4d& a,const Vec4d& b,const Vec4d& c)
{
    return Vec4d(minmod(a.x(),b.x(),c.x()),
                 minmod(a.y(),b.y(),c.y()),
                 minmod(a.z(),b.z(),c.z()),
                 minmod(a.w(),b.w(),c.w()));
}

template<typename T>
FORCE_INLINE T getAngle3D(const typename ScalarUtil<T>::ScalarVec3& a,const typename ScalarUtil<T>::ScalarVec3& b)
{
    T denom=std::max<T>(a.norm()*b.norm(),ScalarUtil<T>::scalar_eps);
    return acos(max<T>(min<T>(a.dot(b)/denom,(T)(1.0f-ScalarUtil<T>::scalar_eps)),(T)(-1.0f+ScalarUtil<T>::scalar_eps)));
}
template<typename T>
FORCE_INLINE T getAngle2D(const typename ScalarUtil<T>::ScalarVec2& a,const typename ScalarUtil<T>::ScalarVec2& b)
{
    T denom=std::max<T>(a.norm()*b.norm(),ScalarUtil<T>::scalar_eps);
    return acos(max<T>(min<T>(a.dot(b)/denom,(T)(1.0f-ScalarUtil<T>::scalar_eps)),(T)(-1.0f+ScalarUtil<T>::scalar_eps)));
}

template<typename T>
FORCE_INLINE typename ScalarUtil<T>::ScalarMat3 makeRotation(const typename ScalarUtil<T>::ScalarVec3& rotation)
{
    typename ScalarUtil<T>::ScalarMat3 ret;//=(typename ScalarUtil<T>::ScalarMat3)::Identity();

    typename ScalarUtil<T>::ScalarVec3 negRot=rotation*-1.0f;
    const T cr = cos( negRot.x() );
    const T sr = sin( negRot.x() );
    const T cp = cos( negRot.y() );
    const T sp = sin( negRot.y() );
    const T cy = cos( negRot.z() );
    const T sy = sin( negRot.z() );

    ret(0,0) = ( cp*cy );
    ret(0,1) = ( cp*sy );
    ret(0,2) = ( -sp );

    const T srsp = sr*sp;
    const T crsp = cr*sp;

    ret(1,0) = ( srsp*cy-cr*sy );
    ret(1,1) = ( srsp*sy+cr*cy );
    ret(1,2) = ( sr*cp );

    ret(2,0) = ( crsp*cy+sr*sy );
    ret(2,1) = ( crsp*sy-sr*cy );
    ret(2,2) = ( cr*cp );

    return ret;
}
template<typename T>
FORCE_INLINE typename ScalarUtil<T>::ScalarVec3 getRotation(const typename ScalarUtil<T>::ScalarMat3& m,bool wrap=true)
{
    const typename ScalarUtil<T>::ScalarVec3 scale(m.block(0,0,1,3).norm(),
            m.block(1,0,1,3).norm(),
            m.block(2,0,1,3).norm());
    const typename ScalarUtil<T>::ScalarVec3 invScale(1.0f/scale.x(),
            1.0f/scale.y(),
            1.0f/scale.z());

    T Y = -asin(m(0,2)*invScale.x());
    const T C = cos(Y);

    T rotx, roty, X, Z;

    if (std::abs(C) >= ScalarUtil<T>::scalar_eps) {
        const T invC = 1.0f/C;
        rotx = m(2,2) * invC * invScale.z();
        roty = m(1,2) * invC * invScale.y();
        X = atan2( roty, rotx ) ;
        rotx = m(0,0) * invC * invScale.x();
        roty = m(0,1) * invC * invScale.y();
        Z = atan2( roty, rotx ) ;
    } else {
        X = 0.0f;
        rotx = m(1,1) * invScale.y();
        roty = -m(1,0) * invScale.y();
        Z = atan2( roty, rotx ) ;
    }

    X*=-1.0f;
    Y*=-1.0f;
    Z*=-1.0f;

    if(wrap) {
        // fix values that get below zero
        // before it would set (!) values to 360
        // that were above 360:
        if (X < 0.0f) X += ((T)M_PI)*2.0f;
        if (Y < 0.0f) Y += ((T)M_PI)*2.0f;
        if (Z < 0.0f) Z += ((T)M_PI)*2.0f;
    }

    return typename ScalarUtil<T>::ScalarVec3(X,Y,Z);
}

template<typename T>
FORCE_INLINE bool getProjectionMatrixFrustum(const typename ScalarUtil<T>::ScalarMat4& prj,
        T& left,T& right,T& bottom,T& top,T& zNear,T& zFar)
{
    if (prj(3,0)!=0.0f || prj(3,1)!=0.0f || prj(3,2)!=-1.0f || prj(3,3)!=0.0f)
        return false;

    // note: near and far must be used inside this method instead of zNear and zFar
    // because zNear and zFar are references and they may point to the same variable.
    T temp_near = prj(2,3) / (prj(2,2)-1.0f);
    T temp_far = prj(2,3) / (1.0f+prj(2,2));

    left = temp_near * (prj(0,2)-1.0f) / prj(0,0);
    right = temp_near * (1.0f+prj(0,2)) / prj(0,0);

    top = temp_near * (1.0f+prj(1,2)) / prj(1,1);
    bottom = temp_near * (prj(1,2)-1.0f) / prj(1,1);

    zNear = temp_near;
    zFar = temp_far;

    return true;
}
template<typename T>
FORCE_INLINE void getViewMatrixFrame(const typename ScalarUtil<T>::ScalarMat4& m,
                                     typename ScalarUtil<T>::ScalarVec3& c,
                                     typename ScalarUtil<T>::ScalarVec3& x,
                                     typename ScalarUtil<T>::ScalarVec3& y,
                                     typename ScalarUtil<T>::ScalarVec3& z)
{
    typename ScalarUtil<T>::ScalarMat4 mInv=m.inverse();

    typename ScalarUtil<T>::ScalarVec4 cH=mInv*typename ScalarUtil<T>::ScalarVec4(0.0f,0.0f,0.0f,1.0f);
    c=(typename ScalarUtil<T>::ScalarVec3(cH.x(),cH.y(),cH.z())*(1.0f/cH.w()));

    typename ScalarUtil<T>::ScalarVec4 xH=mInv*typename ScalarUtil<T>::ScalarVec4(1.0f,0.0f,0.0f,1.0f);
    x=(typename ScalarUtil<T>::ScalarVec3(xH.x(),xH.y(),xH.z())*(1.0f/xH.w())-c).normalized();

    typename ScalarUtil<T>::ScalarVec4 yH=mInv*typename ScalarUtil<T>::ScalarVec4(0.0f,1.0f,0.0f,1.0f);
    y=(typename ScalarUtil<T>::ScalarVec3(yH.x(),yH.y(),yH.z())*(1.0f/yH.w())-c).normalized();

    typename ScalarUtil<T>::ScalarVec4 zH=mInv*typename ScalarUtil<T>::ScalarVec4(0.0f,0.0f,-1.0f,1.0f);
    z=(typename ScalarUtil<T>::ScalarVec3(zH.x(),zH.y(),zH.z())*(1.0f/zH.w())-c).normalized();
}
template<typename T>
FORCE_INLINE bool getViewMatrixFrame(const typename ScalarUtil<T>::ScalarMat4& m,
                                     const typename ScalarUtil<T>::ScalarMat4& prj,
                                     typename ScalarUtil<T>::ScalarVec3& c,
                                     typename ScalarUtil<T>::ScalarVec3& x,
                                     typename ScalarUtil<T>::ScalarVec3& y,
                                     typename ScalarUtil<T>::ScalarVec3& z)
{
    T left,right,bottom,top,zNear,zFar;
    if(!getProjectionMatrixFrustum<T>(prj,left,right,bottom,top,zNear,zFar))
        return false;

    typename ScalarUtil<T>::ScalarMat4 mInv=m.inverse();
    typename ScalarUtil<T>::ScalarMat3 mRotInv=(m.block(0,0,3,3)).inverse();

    typename ScalarUtil<T>::ScalarVec4 cH=mInv*typename ScalarUtil<T>::ScalarVec4(0.0f,0.0f,0.0f,1.0f);
    c=(typename ScalarUtil<T>::ScalarVec3(cH.x(),cH.y(),cH.z())*(1.0f/cH.w()));

    x=mRotInv*(typename ScalarUtil<T>::ScalarVec3(1.0f,0.0f,0.0f));
    x=x.normalized()*right;

    y=mRotInv*(typename ScalarUtil<T>::ScalarVec3(0.0f,1.0f,0.0f));
    y=y.normalized()*top;

    z=mRotInv*(typename ScalarUtil<T>::ScalarVec3(0.0f,0.0f,-1.0f));
    z=z.normalized()*zNear;

    return true;
}
template<typename T>
FORCE_INLINE typename ScalarUtil<T>::ScalarVec3
transformHomo(const typename ScalarUtil<T>::ScalarMat4& m,
              const typename ScalarUtil<T>::ScalarVec3& p)
{
    const T B=(m(3,0)*p.x()+m(3,1)*p.y()+m(3,2)*p.z()+m(3,3));
    const typename ScalarUtil<T>::ScalarVec3 A(m(0,0)*p.x()+m(0,1)*p.y()+m(0,2)*p.z()+m(0,3),
            m(1,0)*p.x()+m(1,1)*p.y()+m(1,2)*p.z()+m(1,3),
            m(2,0)*p.x()+m(2,1)*p.y()+m(2,2)*p.z()+m(2,3));
    return A/B;
}
template<typename T>
FORCE_INLINE typename ScalarUtil<T>::ScalarMat4
getOrtho2D(const typename ScalarUtil<T>::ScalarVec3& minC,
           const typename ScalarUtil<T>::ScalarVec3& maxC)
{
    typename ScalarUtil<T>::ScalarMat4 mt=typename ScalarUtil<T>::ScalarMat4::Zero();
    mt(0,0)=2.0f/(maxC.x()-minC.x());
    mt(0,3)=(minC.x()+maxC.x())/(minC.x()-maxC.x());

    mt(1,1)=2.0f/(maxC.y()-minC.y());
    mt(1,3)=(minC.y()+maxC.y())/(minC.y()-maxC.y());

    mt(2,2)=1.0f;
    mt(3,3)=1.0f;
    return mt;
}

//linear interpolation
template <typename T,typename TF>
FORCE_INLINE T interp1D(const T& v0,const T& v1,const TF& px)
{
    return v0*(1.0f-px)+v1*px;
}
template <typename T,typename TF>
FORCE_INLINE T interp2D(const T& v0,const T& v1,const T& v2,const T& v3,const TF& px,const TF& py)
{
    return interp1D(interp1D(v0,v1,px),
                    interp1D(v2,v3,px),py);
}
template <typename T,typename TF>
FORCE_INLINE T interp3D(const T& v0,const T& v1,const T& v2,const T& v3,
                        const T& v4,const T& v5,const T& v6,const T& v7,
                        const TF& px,const TF& py,const TF& pz)
{
    return interp1D(interp2D(v0,v1,v2,v3,px,py),
                    interp2D(v4,v5,v6,v7,px,py),pz);
}

//inter interpolation stencil
template <typename TF>
FORCE_INLINE sizeType stencil1D(TF* coefs,const TF& px)
{
    coefs[0]=1.0f-px;
    coefs[1]=px;
    return 2;
}
template <typename TF>
FORCE_INLINE sizeType stencil2D(TF* coefs,const TF& px,const TF& py)
{
    coefs[0]=(1.0f-px)*(1.0f-py);
    coefs[1]=px*(1.0f-py);
    coefs[2]=(1.0f-px)*py;
    coefs[3]=px*py;
    return 4;
}
template <typename TF>
FORCE_INLINE sizeType stencil3D(TF* coefs,const TF& px,const TF& py,const TF& pz)
{
    coefs[0]=(1.0f-px)*(1.0f-py)*(1.0f-pz);
    coefs[1]=px*(1.0f-py)*(1.0f-pz);
    coefs[2]=(1.0f-px)*py*(1.0f-pz);
    coefs[3]=px*py*(1.0f-pz);
    coefs[4]=(1.0f-px)*(1.0f-py)*pz;
    coefs[5]=px*(1.0f-py)*pz;
    coefs[6]=(1.0f-px)*py*pz;
    coefs[7]=px*py*pz;
    return 8;
}

//linear interpolation with minmax
template <typename T,typename TF>
FORCE_INLINE T interp1D(const T& v0,const T& v1,const TF& px,T& minV,T& maxV)
{
    minV=compMin(v0,v1);
    maxV=compMax(v0,v1);
    return v0*(1.0f-px)+v1*px;
}
template <typename T,typename TF>
FORCE_INLINE T interp2D(const T& v0,const T& v1,const T& v2,const T& v3,const TF& px,const TF& py,T& minV,T& maxV)
{
    T minV2,maxV2;
    T ret=interp1D(interp1D(v0,v1,px,minV,maxV),
                   interp1D(v2,v3,px,minV2,maxV2),py);
    minV=compMin(minV,minV2);
    maxV=compMax(maxV,maxV2);
    return ret;
}
template <typename T,typename TF>
FORCE_INLINE T interp3D(const T& v0,const T& v1,const T& v2,const T& v3,
                        const T& v4,const T& v5,const T& v6,const T& v7,
                        const TF& px,const TF& py,const TF& pz,T& minV,T& maxV)
{
    T minV2,maxV2;
    T ret=interp1D(interp2D(v0,v1,v2,v3,px,py,minV,maxV),
                   interp2D(v4,v5,v6,v7,px,py,minV2,maxV2),pz);
    minV=compMin(minV,minV2);
    maxV=compMax(maxV,maxV2);
    return ret;
}

#ifdef _MSC_VER
#define STRINGIFY(X) X
#define PRAGMA __pragma
#else
#define STRINGIFY(X) #X
#define PRAGMA _Pragma
#endif

#define OMP_PARALLEL_FOR_ PRAGMA(STRINGIFY(omp parallel for num_threads(OmpSettings::getOmpSettings().nrThreads()) schedule(static,OmpSettings::getOmpSettings().szChunk())))
#define OMP_PARALLEL_FOR_I(...) PRAGMA(STRINGIFY(omp parallel for num_threads(OmpSettings::getOmpSettings().nrThreads()) schedule(static,OmpSettings::getOmpSettings().szChunk()) __VA_ARGS__))
#define OMP_PARALLEL_FOR_X(X) PRAGMA(STRINGIFY(omp parallel for num_threads(X) schedule(static,OmpSettings::getOmpSettings().szChunk())))
#define OMP_ADD(...) reduction(+:##__VA_ARGS__)
#define OMP_PRI(...) private(__VA_ARGS__)
#define OMP_FPRI(...) firstprivate(__VA_ARGS__)
#define OMP_ATOMIC_	PRAGMA(STRINGIFY(omp atomic))
#define OMP_ATOMIC_CAPTURE_	PRAGMA(STRINGIFY(omp atomic capture))
#define OMP_CRITICAL_ PRAGMA(STRINGIFY(omp critical))
#define OMP_FLUSH_(X) PRAGMA(STRINGIFY(omp flush(X)))
struct OmpSettings {
public:
    static const OmpSettings& getOmpSettings() {return _ompSettings;}
    int nrThreads() const {return _nrThreads;}
    int szChunk() const {return _szChunk;}
private:
    OmpSettings():_nrThreads(std::max<int>(omp_get_num_procs(),2)*3/4),_szChunk(4) {}
    const int _nrThreads;
    const int _szChunk;
    static OmpSettings _ompSettings;
};

// #define TIME_PERFORMANCE
#ifndef TIME_PERFORMANCE
#define TBEG(STR) {
#define TEND }

#define TBEGT {
#define TENDT(STR,tree) }
#else
#define TBEG(STR) {std::cout << STR << ": ";boost::timer::cpu_timer __t;__t.start();
#define TEND __t.stop();std::cout << __t.format() << std::endl;}

#define TBEGT {boost::timer::cpu_timer __t;__t.start();
#define TENDT(STR,tree) __t.stop();tree.put<scalarD>(STR,__t.elapsed().wall/1000000000.0L);}
#endif

PRJ_END

#endif
