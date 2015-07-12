#ifndef IMPLICIT_FUNC_INTERFACE_H
#define IMPLICIT_FUNC_INTERFACE_H

#include "MathBasic.h"

PRJ_BEGIN

template <typename T>
class ImplicitFunc
{
public:
    virtual ~ImplicitFunc() {}
    virtual T operator()(const typename ScalarUtil<T>::ScalarVec3& pos) const =0;
};

template <typename T>
class VelFunc
{
public:
    virtual ~VelFunc() {}
    virtual typename ScalarUtil<T>::ScalarVec3 operator()(const typename ScalarUtil<T>::ScalarVec3& pos) const =0;
};

PRJ_END

#endif