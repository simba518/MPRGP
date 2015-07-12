#ifndef CALLBACK_H
#define CALLBACK_H

#include "MathBasic.h"

PRJ_BEGIN

//============================================================================
// Iteration callback
template <typename T,typename KERNEL_TYPE>
struct Callback {
public:
    typedef typename KERNEL_TYPE::Vec Vec;
    virtual ~Callback() {}
    virtual void reset() {
        //INFOV("Start")
    }
    virtual sizeType operator()(const Vec& x,const T& residueNorm,const sizeType& n) {
        INFOV("ITER %d: r=%f",(int)n,residueNorm)
        return 0;
    }
    virtual sizeType operator()(const Vec& x,const Vec& g,const T& fx,const T& xnorm,const T& gnorm,const T& step,const sizeType& k,const sizeType& ls) {
        INFOV("ITER %d: f=%f, x=%f, g=%f, s=%f",(int)k,fx,xnorm,gnorm,step)
        return 0;
    }
};

template <typename T,typename KERNEL_TYPE>
struct CallbackCounter : public Callback<T,KERNEL_TYPE> {
public:
    typedef typename Callback<T,KERNEL_TYPE>::Vec Vec;
    CallbackCounter():_n(0) {}
    virtual void reset() {
        INFOV("Start %d",_n++)
    }
    virtual sizeType operator()(const Vec& x,const T& residueNorm,const sizeType& n) {
        return 0;
    }
    virtual sizeType operator()(const Vec& x,const Vec& g,const T& fx,const T& xnorm,const T& gnorm,const T& step,const sizeType& k,const sizeType& ls) {
        return 0;
    }
    sizeType _n;
};

template <typename T,typename KERNEL_TYPE>
struct ConvergencyCounter : public Callback<T,KERNEL_TYPE> {
public:
    typedef typename Callback<T,KERNEL_TYPE>::Vec Vec;
	ConvergencyCounter(const std::string& path):_os(path.c_str(),std::ios::out){}
	~ConvergencyCounter()
	{
		for(sizeType i=0;i<(sizeType)_conv.size();i++)
			_os << i+1 << " " << _conv[i].first/_conv[i].second << std::endl;
	}
	virtual void reset(){_ind=0;}
    virtual sizeType operator()(const Vec& x,const T& residueNorm,const sizeType& n)
	{
		if(_ind == 0){
			_initRes=residueNorm;
		}else{
			if((sizeType)_conv.size() < _ind)
				_conv.push_back(make_pair(residueNorm/_initRes,1));
			else{
				_conv[_ind-1].first+=residueNorm/_initRes;
				_conv[_ind-1].second+=1;
			}
		}
		_ind++;
		return 0;
	}
protected:
	std::ofstream _os;
	vector<std::pair<T,sizeType> > _conv;
	sizeType _ind;
	T _initRes;
};

PRJ_END

#endif
