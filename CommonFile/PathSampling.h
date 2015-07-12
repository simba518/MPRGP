#ifndef PATH_SAMPLING_H
#define PATH_SAMPLING_H

#include "MathBasic.h"

PRJ_BEGIN

template <typename VEC>
class PathSampling
{
    static const int dim=VEC::RowsAtCompileTime;
    typedef typename VEC::Scalar T;
public:
    PathSampling(const std::vector<VEC,Eigen::aligned_allocator<VEC> >& path)
        :_path(path),_frac(0.0f) {}
    virtual ~PathSampling() {}
    bool isEnded() const {
        return _frac >= (T)(_path.size()-1);
    }
    VEC getPoint() const {
        return interp(_path);
    }
    VEC interp(const std::vector<VEC,Eigen::aligned_allocator<VEC> >& data) const {
        if(isEnded())
            return data.back();
        else {
            const sizeType base=(sizeType)_frac;
            const T frac=_frac-(T)base;
            return data[base]*(1.0f-frac)+data[base+1]*frac;
        }
    }
    bool advance(const T& len) {
        T totalDist=0.0f;
        while(_frac < (scalar)_path.size()-1.0f && totalDist < len) {
            const sizeType base=(sizeType)_frac;
            const T frac=_frac-base;
            const T currDist=(_path[base+1]-_path[base]).norm()*(1.0f-frac);

            if(totalDist+currDist < len) {
                _frac=(scalar)base+1.0f;
                totalDist+=currDist;
            } else {
                _frac+=(1.0f-frac)*(len-totalDist)/(currDist+EPS);
                totalDist=len;
            }
        }
        return isEnded();
    }
protected:
    const std::vector<VEC,Eigen::aligned_allocator<VEC> >& _path;
    T _frac;
};

PRJ_END

#endif
