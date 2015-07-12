#ifndef _HARMONICOSCILLATOR_H_
#define _HARMONICOSCILLATOR_H_

#include <math.h>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <assertext.h>

namespace UTILITY{
  
  /**
   * @class HarmonicOscillator analysis solution for the motion equation:
   * \ddot{z} + d\dot{z} + \lambda z = 0.
   * where usually we set d = \alpha_k\lambda + \alpha_m.
   * @see Interactive Shape Interpolation through Controllable Dynamic Deformation.
   */
  template <class T>
  class HarmonicOscillator{
	
  public:
	HarmonicOscillator(const T lambda,const T d,const T z0,const T v0):
	  _lambda(lambda),_d(d),_z0(z0),_v0(v0){
	  precompute();
	}

	HarmonicOscillator(const T lambda,const T alpha_k,const T alpha_m,const T z0,const T v0):
	  _lambda(lambda),_d(alpha_k*lambda+alpha_m),_z0(z0),_v0(v0){
	  assert_ge(alpha_k,0);
	  assert_ge(alpha_m,0);
	  precompute();
	}

	void precompute(){
	  assert_gt(_lambda,0);
	  assert_ge(_d,0);
	  _alpha = _d/2;
	  assert_gt((_lambda-_alpha*_alpha), 0);
	  _w = sqrt(_lambda - _alpha*_alpha);
	  _P = _z0;
	  _Q = (_v0+_alpha*_z0)/_w;
	}

	T operator() (const T t)const{
	  assert_ge(t,0);
	  return (_P*cos(_w*t) + _Q*sin(_w*t))*exp(-1.0f*_alpha*t);
	}

	template <typename VECTOR>
	const VECTOR &generateSequence(const T t0,const T dt,const int len,VECTOR &q)const{

	  assert_ge(t0,0);
	  assert_gt(dt,0);
	  assert_ge(len,0);
	  q.resize(len);
	  for (int i = 0; i < len; ++i){
		q[i] = (*this)(t0+((double)i)*dt);
	  }
	  return q;
	}

	template <typename VECTOR> 
	VECTOR generateSequence(const T t0,const T dt,const int len)const{
	  VECTOR q;
	  generateSequence(t0,dt,len,q);
	  return q;
	}
	
  private:
	T _lambda;
	T _d;
	T _z0;
	T _v0;

	T _P;
	T _Q;
	T _w;
	T _alpha;
  };
  
  typedef boost::shared_ptr<HarmonicOscillator<double> > pHarmonicOscillatorD;
  typedef boost::shared_ptr<HarmonicOscillator<float> > pHarmonicOscillatorF;

  template <class T>
  class HarmonicOscillatorSet{
	
  public:
	template <typename VECTOR>
	HarmonicOscillatorSet(const VECTOR&lambda,const T d,const VECTOR&z0,const VECTOR&v0){

	  assert_eq(lambda.size(),z0.size());
	  assert_eq(lambda.size(),v0.size());
	  _oscillators.clear();
	  for (size_t i = 0; i < lambda.size(); ++i){
		_oscillators.push_back(HarmonicOscillator<T>(lambda[i],d,z0[i],v0[i]));
		_oscillators[i].precompute();
	  }
	}

	template <typename VECTOR>
	HarmonicOscillatorSet(const VECTOR&lambda,const VECTOR d,const VECTOR&z0,const VECTOR&v0){

	  assert_eq(lambda.size(),z0.size());
	  assert_eq(lambda.size(),v0.size());
	  assert_eq(lambda.size(),d.size());
	  _oscillators.clear();
	  for (size_t i = 0; i < lambda.size(); ++i){
		_oscillators.push_back(HarmonicOscillator<T>(lambda[i],d[i],z0[i],v0[i]));
		_oscillators[i].precompute();
	  }
	}

	template <typename VECTOR>
	HarmonicOscillatorSet(const VECTOR &lambda,const T alpha_k,const T alpha_m,const VECTOR &z0,const VECTOR &v0){
	  assert_eq(lambda.size(),z0.size());
	  assert_eq(lambda.size(),v0.size());
	  _oscillators.clear();
	  for (int i = 0; i < (int)lambda.size(); ++i){
		_oscillators.push_back(HarmonicOscillator<T>(lambda[i],alpha_k,alpha_m,z0[i],v0[i]));
		_oscillators[i].precompute();
	  }
	}

	template <typename VECTOR>
	VECTOR at(const T t)const{
	  assert_ge(t,0);
	  VECTOR z(_oscillators.size());
	  for (int i = 0; i < z.size(); ++i)
		z[i] = _oscillators[i](t);
	  return z;
	}

	template <typename MATRIX>
	MATRIX generateSequence(const T t0,const T dt,const int len)const{

	  MATRIX Z(_oscillators.size(),len);
	  std::vector<T> q;
	  for (size_t i = 0; i < _oscillators.size(); ++i){
		_oscillators[i].generateSequence(t0,dt,len,q);
		for (size_t j = 0; j < q.size(); ++j){
		  Z(i,j) = q[j];
		}
	  }
	  return Z;
	}
	
  private:
	std::vector<HarmonicOscillator<T> > _oscillators;
  };
  
}//end of namespace

#endif /*_HARMONICOSCILLATOR_H_*/
