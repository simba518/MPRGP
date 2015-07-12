#ifndef DATA_CONVERT_H
#define DATA_CONVERT_H

#include "MathBasic.h"

PRJ_BEGIN

template <typename T1,typename T2> 
FORCE_INLINE void transfer(const Eigen::Matrix<T1,2,1>& from,Eigen::Matrix<T2,2,1>& to)
{
    to[0]=(T2)from[0];
    to[1]=(T2)from[1];
}
template <typename VEC,typename T2> 
FORCE_INLINE void transfer(const Eigen::Block<VEC,2,1>& from,Eigen::Matrix<T2,2,1>& to)
{
    to[0]=(T2)from[0];
    to[1]=(T2)from[1];
}
template <typename T1,typename VEC> 
FORCE_INLINE void transfer(const Eigen::Matrix<T1,2,1>& from,Eigen::Block<VEC,2,1>& to)
{
    typedef typename VEC::Scalar TO;
    to[0]=(TO)from[0];
    to[1]=(TO)from[1];
}

template <typename T1,typename T2> 
FORCE_INLINE void transfer(const Eigen::Matrix<T1,3,1>& from,Eigen::Matrix<T2,3,1>& to)
{
    to[0]=(T2)from[0];
    to[1]=(T2)from[1];
    to[2]=(T2)from[2];
}
template <typename VEC,typename T2> 
FORCE_INLINE void transfer(const Eigen::Block<VEC,3,1>& from,Eigen::Matrix<T2,3,1>& to)
{
    to[0]=(T2)from[0];
    to[1]=(T2)from[1];
    to[2]=(T2)from[2];
}
template <typename T1,typename VEC> 
FORCE_INLINE void transfer(const Eigen::Matrix<T1,3,1>& from,Eigen::Block<VEC,3,1>& to)
{
    typedef typename VEC::Scalar TO;
    to[0]=(TO)from[0];
    to[1]=(TO)from[1];
    to[2]=(TO)from[2];
}

template <typename T1,typename T2> 
FORCE_INLINE void transfer(const Eigen::Matrix<T1,4,1>& from,Eigen::Matrix<T2,4,1>& to)
{
    to[0]=(T2)from[0];
    to[1]=(T2)from[1];
    to[2]=(T2)from[2];
    to[3]=(T2)from[3];
}
template <typename VEC,typename T2> 
FORCE_INLINE void transfer(const Eigen::Block<VEC,4,1>& from,Eigen::Matrix<T2,4,1>& to)
{
    to[0]=(T2)from[0];
    to[1]=(T2)from[1];
    to[2]=(T2)from[2];
    to[3]=(T2)from[3];
}
template <typename T1,typename VEC> 
FORCE_INLINE void transfer(const Eigen::Matrix<T1,4,1>& from,Eigen::Block<VEC,4,1>& to)
{
    typedef typename VEC::Scalar TO;
    to[0]=(TO)from[0];
    to[1]=(TO)from[1];
    to[2]=(TO)from[2];
    to[3]=(TO)from[3];
}

template <typename T1,typename T2> 
FORCE_INLINE void transfer(const Eigen::Quaternion<T1>& from,Eigen::Quaternion<T2>& to)
{
    to.w()=(T2)from.w();
    to.x()=(T2)from.x();
    to.y()=(T2)from.y();
    to.z()=(T2)from.z();
}

template <typename VEC,typename T> 
FORCE_INLINE void transfer(const VEC& from,Eigen::Quaternion<T>& to)
{
    to.w()=(T)from[0];
    to.x()=(T)from[1];
    to.y()=(T)from[2];
    to.z()=(T)from[3];
}
template <typename VEC,typename T> 
FORCE_INLINE void transfer(const Eigen::Block<VEC,4,1>& from,Eigen::Quaternion<T>& to)
{
    to.w()=(T)from[0];
    to.x()=(T)from[1];
    to.y()=(T)from[2];
    to.z()=(T)from[3];
}

template <typename T,typename VEC> 
FORCE_INLINE void transfer(const Eigen::Quaternion<T>& from,VEC& to)
{
    typedef typename VEC::Scalar TO;
    to[0]=(TO)from.w();
    to[1]=(TO)from.x();
    to[2]=(TO)from.y();
    to[3]=(TO)from.z();
}
template <typename T,typename VEC> 
FORCE_INLINE void transfer(const Eigen::Quaternion<T>& from,Eigen::Block<VEC,4,1>& to)
{
    typedef typename VEC::Scalar TO;
    to[0]=(TO)from.w();
    to[1]=(TO)from.x();
    to[2]=(TO)from.y();
    to[3]=(TO)from.z();
}

template <typename T1,typename T2> 
FORCE_INLINE void transfer(const Eigen::Matrix<T1,3,3>& from,Eigen::Matrix<T2,3,3>& to)
{
    to(0,0)=(T2)from(0,0);
    to(0,1)=(T2)from(0,1);
    to(0,2)=(T2)from(0,2);
    //to(0,3)=(T2)from(0,3);

    to(1,0)=(T2)from(1,0);
    to(1,1)=(T2)from(1,1);
    to(1,2)=(T2)from(1,2);
    //to(1,3)=(T2)from(1,3);

    to(2,0)=(T2)from(2,0);
    to(2,1)=(T2)from(2,1);
    to(2,2)=(T2)from(2,2);
    //to(2,3)=(T2)from(2,3);
}
template <typename T1,typename T2> 
FORCE_INLINE void transfer(const Eigen::Matrix<T1,4,4>& from,Eigen::Matrix<T2,4,4>& to)
{
    to(0,0)=(T2)from(0,0);
    to(0,1)=(T2)from(0,1);
    to(0,2)=(T2)from(0,2);
    to(0,3)=(T2)from(0,3);

    to(1,0)=(T2)from(1,0);
    to(1,1)=(T2)from(1,1);
    to(1,2)=(T2)from(1,2);
    to(1,3)=(T2)from(1,3);

    to(2,0)=(T2)from(2,0);
    to(2,1)=(T2)from(2,1);
    to(2,2)=(T2)from(2,2);
    to(2,3)=(T2)from(2,3);

    to(3,0)=(T2)from(3,0);
    to(3,1)=(T2)from(3,1);
    to(3,2)=(T2)from(3,2);
    to(3,3)=(T2)from(3,3);
}
template <typename T1,typename T2> 
FORCE_INLINE void transfer(const BBox<T1>& from,BBox<T2>& to)
{
    transfer(from._minC,to._minC);
    transfer(from._maxC,to._maxC);
}
template <typename FROM,typename TO> FORCE_INLINE TO transferOut(const FROM& from){
    TO ret;transfer(from,ret);return ret;
}

PRJ_END

#endif