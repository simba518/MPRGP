#ifndef MINIMIZER_H
#define MINIMIZER_H

#include "MathBasic.h"
#include "Callback.h"
#include "Objective.h"
#include "IO.h"
#include <vector>

PRJ_BEGIN

//Just a port of the code by:
//Copyright (c) 1990, Jorge Nocedal
//Copyright (c) 2007-2010 Naoaki Okazaki
enum STATUS_CODE {
    /** L-BFGS reaches convergence. */
    SUCCESS = 0,
    CONVERGENCE = 0,
    STOP,
    /** The initial variables already minimize the objective function. */
    ALREADY_MINIMIZED,

    /** Unknown error. */
    ERR_UNKNOWNERROR = -1024,
    /** Logic error. */
    ERR_LOGICERROR,
    /** The minimization process has been canceled. */
    ERR_CANCELED,
    /** Invalid number of variables specified. */
    ERR_INVALID_N,
    /** Invalid number of variables (for SSE) specified. */
    ERR_INVALID_N_SSE,
    /** The array x must be aligned to 16 (for SSE). */
    ERR_INVALID_X_SSE,
    /** Invalid parameter lbfgs_parameter_t::epsilon specified. */
    ERR_INVALID_EPSILON,
    /** Invalid parameter lbfgs_parameter_t::past specified. */
    ERR_INVALID_TESTPERIOD,
    /** Invalid parameter lbfgs_parameter_t::delta specified. */
    ERR_INVALID_DELTA,
    /** Invalid parameter lbfgs_parameter_t::linesearch specified. */
    ERR_INVALID_LINESEARCH,
    /** Invalid parameter lbfgs_parameter_t::max_step specified. */
    ERR_INVALID_MINSTEP,
    /** Invalid parameter lbfgs_parameter_t::max_step specified. */
    ERR_INVALID_MAXSTEP,
    /** Invalid parameter lbfgs_parameter_t::ftol specified. */
    ERR_INVALID_FTOL,
    /** Invalid parameter lbfgs_parameter_t::wolfe specified. */
    ERR_INVALID_WOLFE,
    /** Invalid parameter lbfgs_parameter_t::gtol specified. */
    ERR_INVALID_GTOL,
    /** Invalid parameter lbfgs_parameter_t::xtol specified. */
    ERR_INVALID_XTOL,
    /** Invalid parameter lbfgs_parameter_t::max_linesearch specified. */
    ERR_INVALID_MAXLINESEARCH,
    /** Invalid parameter lbfgs_parameter_t::orthantwise_c specified. */
    ERR_INVALID_ORTHANTWISE,
    /** Invalid parameter lbfgs_parameter_t::orthantwise_start specified. */
    ERR_INVALID_ORTHANTWISE_START,
    /** Invalid parameter lbfgs_parameter_t::orthantwise_end specified. */
    ERR_INVALID_ORTHANTWISE_END,
    /** The line-search step went out of the interval of uncertainty. */
    ERR_OUTOFINTERVAL,
    /** A logic error occurred; alternatively, the interval of uncertainty
    	became too small. */
    ERR_INCORRECT_TMINMAX,
    /** A rounding error occurred; alternatively, no line-search step
    	satisfies the sufficient decrease and curvature conditions. */
    ERR_ROUNDING_ERROR,
    /** The line-search step became smaller than lbfgs_parameter_t::min_step. */
    ERR_MINIMUMSTEP,
    /** The line-search step became larger than lbfgs_parameter_t::max_step. */
    ERR_MAXIMUMSTEP,
    /** The line-search routine reaches the maximum number of evaluations. */
    ERR_MAXIMUMLINESEARCH,
    /** The algorithm routine reaches the maximum number of iterations. */
    ERR_MAXIMUMITERATION,
    /** Relative width of the interval of uncertainty is at most
    	lbfgs_parameter_t::xtol. */
    ERR_WIDTHTOOSMALL,
    /** A logic error (negative line-search step) occurred. */
    ERR_INVALIDPARAMETERS,
    /** The current search direction increases the objective function value. */
    ERR_INCREASEGRADIENT,
    /** User want to cancel. */
    ERR_USER_ASKED,
};
enum LINE_SEARCHER {
    LINESEARCH_DEFAULT = 0,
    LINESEARCH_MORETHUENTE = 0,
    LINESEARCH_BACKTRACKING_ARMIJO = 1,
    LINESEARCH_BACKTRACKING = 2,
    LINESEARCH_BACKTRACKING_WOLFE = 2,
    LINESEARCH_BACKTRACKING_STRONG_WOLFE = 3,
};
enum PROBLEM_PROFILE_CODE {
    PROFILE_LINE_SEARCH=1,
};

template < typename T, typename KERNEL_TYPE=Kernel<T> >
class LineSearcher
{
public:
    typedef typename KERNEL_TYPE::Vec Vec;
    LineSearcher() {
        _maxLineSearch=40;
        _minStep=1e-20f;
        _maxStep=1e20f;
        _ftol=1e-4f;
        _wolfe=0.9f;
        _gtol=0.9f;
        _xtol=1.0e-16f;
    }
    virtual sizeType check(LINE_SEARCHER searcher) const {
        if(_minStep < 0.0f)
            return ERR_INVALID_MINSTEP;
        if(_maxStep < _minStep)
            return ERR_INVALID_MAXSTEP;
        if(_ftol < 0.0f)
            return ERR_INVALID_FTOL;
        if(searcher == LINESEARCH_BACKTRACKING_WOLFE ||
                searcher == LINESEARCH_BACKTRACKING_STRONG_WOLFE) {
            if(_wolfe <= _ftol || 1.0f <= _wolfe)
                return ERR_INVALID_WOLFE;
        }
        if(_gtol < 0.0f)
            return ERR_INVALID_GTOL;
        if(_xtol < 0.0f)
            return ERR_INVALID_XTOL;
        if(_maxLineSearch <= 0)
            return ERR_INVALID_MAXLINESEARCH;
        return SUCCESS;
    }
    virtual sizeType operator()(Vec& x,T& f,Vec& g,Vec& d,
                                T& stp,const Vec& XP,const Vec& gp,
                                Vec& wp,Objective<T,KERNEL_TYPE>& func,
                                LINE_SEARCHER searcher) {
        ASSERT_MSG(false,"Abstract Class: Use Subclass")
        return 0;
    }
public:
    const sizeType& maxLineSearch() const {
        return _maxLineSearch;
    }
    sizeType& maxLineSearch() {
        return _maxLineSearch;
    }
    const T& minStepSize() const {
        return _minStep;
    }
    T& minStepSize() {
        return _minStep;
    }
    const T& maxStepSize() const {
        return _maxStep;
    }
    T& maxStepSize() {
        return _maxStep;
    }
    const T& lineSearchAccuracyFunc() const {
        return _ftol;
    }
    T& lineSearchAccuracyFunc() {
        return _ftol;
    }
    const T& wolfeLineSearchCondition() const {
        return _wolfe;
    }
    T& wolfeLineSearchCondition() {
        return _wolfe;
    }
    const T& lineSearchAccuracyGrad() const {
        return _gtol;
    }
    T& lineSearchAccuracyGrad() {
        return _gtol;
    }
    const T& lineSearchAccuracyX() const {
        return _xtol;
    }
    T& lineSearchAccuracyX() {
        return _xtol;
    }
protected:
    sizeType _maxLineSearch;
    T _minStep;
    T _maxStep;
    T _ftol;
    T _wolfe;
    T _gtol;
    T _xtol;
};

template < typename T, typename KERNEL_TYPE=Kernel<T> >
class LineSearcherBackTracking : public LineSearcher<T,KERNEL_TYPE>
{
public:
    typedef typename LineSearcher<T,KERNEL_TYPE>::Vec Vec;
    using LineSearcher<T,KERNEL_TYPE>::_maxLineSearch;
    using LineSearcher<T,KERNEL_TYPE>::_minStep;
    using LineSearcher<T,KERNEL_TYPE>::_maxStep;
    using LineSearcher<T,KERNEL_TYPE>::_ftol;
    using LineSearcher<T,KERNEL_TYPE>::_wolfe;
    using LineSearcher<T,KERNEL_TYPE>::_gtol;
    using LineSearcher<T,KERNEL_TYPE>::_xtol;
    sizeType operator()(Vec& x,T& f,Vec& g,Vec& d,
                        T& stp,const Vec& XP,const Vec& gp,
                        Vec& wp,Objective<T,KERNEL_TYPE>& func,
                        LINE_SEARCHER seacher) {
        //const sizeType N=x.size();
        sizeType count=0;
        T width,dg;
        T finit,dginit=0.0f,dgtest;
        const T dec=0.5f,inc=2.1f;

        //Check the input parameters for errors
        if(stp<=0.0f)
            return ERR_INVALIDPARAMETERS;

        //Compute the initial gradient in the search direction
        dginit=KERNEL_TYPE::dot(g,d);

        //Make sure that s points to a descent direction
        if(0 < dginit)
            return ERR_INCREASEGRADIENT;

        //The initial value of the objective function
        finit=f;
        dgtest=_ftol*dginit;

        for(;;) {
            KERNEL_TYPE::copy(XP,x);
            KERNEL_TYPE::addScaled(stp,d,x);

            //Evaluate the function and gradient values
            if(func(x,f,g,stp,true) < 0)
                return ERR_USER_ASKED;

            ++count;

            if(f > finit+stp*dgtest) {
                width=dec;
            } else {
                //The sufficient decrease condition (Armijo condition)
                if(seacher == LINESEARCH_BACKTRACKING_ARMIJO)
                    //Exit with the Armijo condition
                    return count;

                //Check the Wolfe condition
                dg=KERNEL_TYPE::dot(g,d);
                if(dg < _wolfe*dginit)
                    width = inc;
                else {
                    if(seacher == LINESEARCH_BACKTRACKING_WOLFE)
                        //Exit with the regular Wolfe condition
                        return count;

                    //Check the strong Wolfe condition
                    if(dg > -_wolfe*dginit)
                        width=dec;
                    else
                        //Exit with the strong Wolfe condition
                        return count;
                }
            }

            if(stp < _minStep)
                //The step is the minimum value
                return ERR_MINIMUMSTEP;
            if (stp > _maxStep)
                //The step is the maximum value
                return ERR_MAXIMUMSTEP;
            if (_maxLineSearch <= count)
                //Maximum number of iteration
                return ERR_MAXIMUMLINESEARCH;

            stp*=width;
        }
    }
};

template < typename T, typename KERNEL_TYPE=Kernel<T> >
class LineSearcherBackTrackingOWLQN : public LineSearcher<T,KERNEL_TYPE>
{
public:
    typedef typename LineSearcher<T,KERNEL_TYPE>::Vec Vec;
    using LineSearcher<T,KERNEL_TYPE>::_maxLineSearch;
    using LineSearcher<T,KERNEL_TYPE>::_minStep;
    using LineSearcher<T,KERNEL_TYPE>::_maxStep;
    using LineSearcher<T,KERNEL_TYPE>::_ftol;
    using LineSearcher<T,KERNEL_TYPE>::_wolfe;
    using LineSearcher<T,KERNEL_TYPE>::_gtol;
    using LineSearcher<T,KERNEL_TYPE>::_xtol;
    sizeType operator()(Vec& x,T& f,Vec& g,Vec& d,
                        T& stp,const Vec& XP,const Vec& gp,
                        Vec& wp,Objective<T,KERNEL_TYPE>& func,
                        LINE_SEARCHER seacher) {
        const sizeType N=x.size();
        sizeType i,count=0;
        T width=0.5f,norm=0.0f;
        T finit=f,dgtest;

        //Check the input parameters for errors
        if(stp <= 0.0f)
            return ERR_INVALIDPARAMETERS;

        //Choose the orthant for the new point
        for(i=0; i < N; ++i)
            wp[i]=(XP[i] == 0.0f) ? -gp[i] : XP[i];

        for(;;) {
            //Update the current point
            KERNEL_TYPE::copy(XP,x);
            KERNEL_TYPE::addScaled(stp,d,x);

            //The current point is projected onto the orthant
            L1Project(x,wp);

            //Evaluate the function and gradient values
            if(func(x,f,g,stp,true) < 0)
                return ERR_USER_ASKED;

            //Compute the L1 norm of the variables and add it to the object value
            norm=L1Norm(x);
            f+=norm*_orthantwiseC;

            ++count;

            dgtest=0.0f;
            for(i=0; i < N; ++i)
                dgtest+=(x[i]-XP[i])*gp[i];

            if(f <= finit+_ftol*dgtest)
                //The sufficient decrease condition
                return count;

            if(stp < _minStep)
                //The step is the minimum value
                return ERR_MINIMUMSTEP;
            if(stp > _maxStep)
                //The step is the maximum value
                return ERR_MAXIMUMSTEP;
            if(_maxLineSearch <= count)
                //Maximum number of iteration
                return ERR_MAXIMUMLINESEARCH;

            stp*=width;
        }
    }
public:
    const T& L1Panelty() const {
        return _orthantwiseC;
    }
    T& L1Panelty() {
        return _orthantwiseC;
    }
    const sizeType& L1StartIndex() const {
        return _orthantwiseStart;
    }
    sizeType& L1StartIndex() {
        return _orthantwiseStart;
    }
    const sizeType& L1EndIndex() const {
        return _orthantwiseEnd;
    }
    sizeType& L1EndIndex() {
        return _orthantwiseEnd;
    }
protected:
    T L1Norm(const Vec& x) const {
        T norm=0.0f;
        for(sizeType i=_orthantwiseStart; i<_orthantwiseEnd; ++i)
            norm+=std::abs(x[i]);
        return norm;
    }
    void L1Project(Vec& x,const Vec& sign) const {
        for (sizeType i=_orthantwiseStart; i < _orthantwiseEnd; ++i)
            if(x[i]*sign[i] <= 0.0f)
                x[i]=0.0f;
    }
protected:
    T _orthantwiseC;
    sizeType _orthantwiseStart;
    sizeType _orthantwiseEnd;
};

template < typename T, typename KERNEL_TYPE=Kernel<T> >
class LineSearcherMoretHuente : public LineSearcher<T,KERNEL_TYPE>
{
public:
    typedef typename LineSearcher<T,KERNEL_TYPE>::Vec Vec;
    using LineSearcher<T,KERNEL_TYPE>::_maxLineSearch;
    using LineSearcher<T,KERNEL_TYPE>::_minStep;
    using LineSearcher<T,KERNEL_TYPE>::_maxStep;
    using LineSearcher<T,KERNEL_TYPE>::_ftol;
    using LineSearcher<T,KERNEL_TYPE>::_wolfe;
    using LineSearcher<T,KERNEL_TYPE>::_gtol;
    using LineSearcher<T,KERNEL_TYPE>::_xtol;
    sizeType operator()(Vec& x,T& f,Vec& g,Vec& d,
                        T& stp,const Vec& XP,const Vec& gp,
                        Vec& wp,Objective<T,KERNEL_TYPE>& func,
                        LINE_SEARCHER seacher) {
        sizeType N=x.size();
        sizeType count=0;
        sizeType brackt,stage1,uinfo=0;
        T dg;
        T stx, fx, dgx;
        T sty, fy, dgy;
        T fxm, dgxm, fym, dgym, fm, dgm;
        T finit, ftest1, dginit, dgtest;
        T width, prevWidth;
        T stmin, stmax;

        //Check the input parameters for errors
        if(stp <= 0.0f)
            return ERR_INVALIDPARAMETERS;

        //Compute the initial gradient in the search direction
        dginit=KERNEL_TYPE::dot(g,d);

        //Make sure that s points to a descent direction
        if(0.0f < dginit)
            return ERR_INCREASEGRADIENT;

        //Initialize local variables
        brackt=0;
        stage1=1;
        finit=f;
        dgtest=_ftol*dginit;
        width=_maxStep-_minStep;
        prevWidth=2.0f*width;

        //The variables stx, fx, dgx contain the values of the step,
        //function, and directional derivative at the best step.
        //The variables sty, fy, dgy contain the value of the step,
        //function, and derivative at the other endpoint of
        //the interval of uncertainty.
        //The variables stp, f, dg contain the values of the step,
        //function, and derivative at the current step.
        stx=sty=0.0f;
        fx=fy=finit;
        dgx=dgy=dginit;

        for(;;) {
            //Set the minimum and maximum steps to correspond to the
            //present interval of uncertainty.
            if (brackt) {
                stmin=std::min<T>(stx,sty);
                stmax=std::max<T>(stx,sty);
            } else {
                stmin=stx;
                stmax=stp+4.0f*(stp-stx);
            }

            //Clip the step in the range of [stpmin, stpmax]
            if(stp < _minStep)
                stp=_minStep;
            if(_maxStep < stp)
                stp=_maxStep;

            //If an unusual termination is to occur then let
            //stp be the lowest point obtained so far
            if((brackt && ((stp <= stmin || stmax <= stp) || _maxLineSearch <= count+1 || uinfo != 0)) ||
                    (brackt && (stmax-stmin <= _xtol*stmax)))
                stp = stx;

            //Compute the current value of x:
            //	x <- x + (*stp) * s
            KERNEL_TYPE::copy(XP,x);
            KERNEL_TYPE::addScaled(stp,d,x);

            //Evaluate the function and gradient values
            if(func(x,f,g,stp,true) < 0)
                return ERR_USER_ASKED;
            dg=KERNEL_TYPE::dot(g,d);

            ftest1=finit+stp*dgtest;
            ++count;

            //Test for errors and convergence
            if(brackt && ((stp <= stmin || stmax <= stp) || uinfo != 0))
                //Rounding errors prevent further progress
                return ERR_ROUNDING_ERROR;
            if(stp == _maxStep && f <= ftest1 && dg <= dgtest)
                //The step is the maximum value
                return ERR_MAXIMUMSTEP;
            if(stp == _minStep && (ftest1 < f || dgtest <= dg))
                //The step is the minimum value
                return ERR_MINIMUMSTEP;
            if(brackt && (stmax - stmin) <= _xtol*stmax)
                //Relative width of the interval of uncertainty is at most xtol
                return ERR_WIDTHTOOSMALL;
            if(_maxLineSearch <= count)
                //Maximum number of iteration
                return ERR_MAXIMUMLINESEARCH;
            if(f <= ftest1 && std::abs(dg) <= _gtol * (-dginit))
                //The sufficient decrease condition and the directional derivative condition hold
                return count;

            //In the first stage we seek a step for which the modified
            //function has a nonpositive value and nonnegative derivative
            if(stage1 && f <= ftest1 && std::min<T>(_ftol, _gtol)*dginit <= dg)
                stage1 = 0;

            //A modified function is used to predict the step only if
            //we have not obtained a step for which the modified
            //function has a nonpositive function value and nonnegative
            //derivative, and if a lower function value has been
            //obtained but the decrease is not sufficient
            if (stage1 && ftest1 < f && f <= fx) {
                //Define the modified function and derivative values
                fm = f - stp * dgtest;
                fxm = fx - stx * dgtest;
                fym = fy - sty * dgtest;
                dgm = dg - dgtest;
                dgxm = dgx - dgtest;
                dgym = dgy - dgtest;

                //Call update_trial_interval() to update the interval of
                //uncertainty and to compute the new step
                uinfo = updateTrialInterval(stx, fxm, dgxm,
                                            sty, fym, dgym,
                                            stp, fm, dgm,
                                            stmin, stmax, brackt);

                //Reset the function and gradient values for f
                fx = fxm + stx * dgtest;
                fy = fym + sty * dgtest;
                dgx = dgxm + dgtest;
                dgy = dgym + dgtest;
            } else {
                //Call update_trial_interval() to update the interval of
                //uncertainty and to compute the new step
                uinfo=updateTrialInterval(stx, fx, dgx,
                                          sty, fy, dgy,
                                          stp, f, dg,
                                          stmin, stmax, brackt);
            }

            //Force a sufficient decrease in the interval of uncertainty
            if (brackt) {
                if (0.66f * prevWidth <= std::abs(sty - stx))
                    stp=stx + 0.5f*(sty - stx);
                prevWidth=width;
                width=std::abs(sty - stx);
            }
        }
        return ERR_LOGICERROR;
    }
protected:
#define USES_MINIMIZER \
    T a, d, gamma, theta, p, q, r, s;

#define CUBIC_MINIMIZER(cm, u, fu, du, v, fv, dv) \
    d = v - u; \
    theta = (fu - fv) * 3.0f / d + du + dv; \
    p = std::abs(theta); \
    q = std::abs(du); \
    r = std::abs(dv); \
    s = std::max<T>( std::max<T>(p, q), r); \
    /* gamma = s*sqrt((theta/s)**2 - (du/s) * (dv/s)) */ \
    a = theta / s; \
    gamma = s * sqrt(a * a - (du / s) * (dv / s)); \
    if (v < u) gamma = -gamma; \
    p = gamma - du + theta; \
    q = gamma - du + gamma + dv; \
    r = p / q; \
    cm = u + r * d;

#define CUBIC_MINIMIZER2(cm, u, fu, du, v, fv, dv, xmin, xmax) \
    d = v - u; \
    theta = (fu - fv) * 3.0f / d + du + dv; \
    p = std::abs(theta); \
    q = std::abs(du); \
    r = std::abs(dv); \
    s = std::max<T>( std::max<T>(p, q), r); \
    /* gamma = s*sqrt((theta/s)**2 - (du/s) * (dv/s)) */ \
    a = theta / s; \
    gamma = s * sqrt(std::max<T>(0, a * a - (du / s) * (dv / s))); \
    if (u < v) gamma = -gamma; \
    p = gamma - dv + theta; \
    q = gamma - dv + gamma + du; \
    r = p / q; \
    if (r < 0.0f && gamma != 0.0f) { \
        cm = v - r * d; \
    } else if (a < 0) { \
        cm = xmax; \
    } else { \
        cm = xmin; \
    }

#define QUARD_MINIMIZER(qm, u, fu, du, v, fv) \
    a = v - u; \
    qm = u + du / ((fu - fv) / a + du) / 2.0f * a;

#define QUARD_MINIMIZER2(qm, u, du, v, dv) \
    a = u - v; \
    qm = v + dv / (dv - du) * a;

#define fsigndiff(x, y) (x * (y / std::abs(y)) < 0.0f)

    sizeType updateTrialInterval(T& x,T& fx,T& dx,
                                 T& y,T& fy,T& dy,
                                 T& t,T& ft,T& dt,
                                 const T tmin,const T tmax,
                                 sizeType& brackt) {
        sizeType bound;
        sizeType dsign=fsigndiff(dt, dx);
        T mc; //minimizer of an interpolated cubic
        T mq; //minimizer of an interpolated quadratic
        T newt;   //new trial value
        USES_MINIMIZER;     //for CUBIC_MINIMIZER and QUARD_MINIMIZER

        //Check the input parameters for errors
        if(brackt) {
            if(t <= std::min<T>(x, y) || std::max<T>(x, y) <= t)
                //The trival value t is out of the interval
                return ERR_OUTOFINTERVAL;
            if (0.0f <= dx*(t - x))
                //The function must decrease from x
                return ERR_INCREASEGRADIENT;
            if(tmax < tmin)
                //Incorrect tmin and tmax specified
                return ERR_INCORRECT_TMINMAX;
        }

        //Trial value selection
        if(fx < ft) {
            //Case 1: a higher function value.
            //The minimum is brackt. If the cubic minimizer is closer
            //to x than the quadratic one, the cubic one is taken, else
            //the average of the minimizers is taken.
            brackt = 1;
            bound = 1;
            CUBIC_MINIMIZER(mc, x, fx, dx, t, ft, dt);
            QUARD_MINIMIZER(mq, x, fx, dx, t, ft);
            if(std::abs(mc - x) < std::abs(mq - x)) {
                newt = mc;
            } else {
                newt = mc + 0.5f * (mq - mc);
            }
        } else if(dsign) {
            //Case 2: a lower function value and derivatives of
            //opposite sign. The minimum is brackt. If the cubic
            //minimizer is closer to x than the quadratic (secant) one,
            //the cubic one is taken, else the quadratic one is taken.
            brackt = 1;
            bound = 0;
            CUBIC_MINIMIZER(mc, x, fx, dx, t, ft, dt);
            QUARD_MINIMIZER2(mq, x, dx, t, dt);
            if(std::abs(mc - t) > std::abs(mq - t)) {
                newt = mc;
            } else {
                newt = mq;
            }
        } else if(std::abs(dt) < std::abs(dx)) {
            //Case 3: a lower function value, derivatives of the
            //same sign, and the magnitude of the derivative decreases.
            //The cubic minimizer is only used if the cubic tends to
            //infinity in the direction of the minimizer or if the minimum
            //of the cubic is beyond t. Otherwise the cubic minimizer is
            //defined to be either tmin or tmax. The quadratic (secant)
            //minimizer is also computed and if the minimum is brackt
            //then the the minimizer closest to x is taken, else the one
            //farthest away is taken.
            bound=1;
            CUBIC_MINIMIZER2(mc, x, fx, dx, t, ft, dt, tmin, tmax);
            QUARD_MINIMIZER2(mq, x, dx, t, dt);
            if(brackt) {
                if(std::abs(t - mc) < std::abs(t - mq)) {
                    newt = mc;
                } else {
                    newt = mq;
                }
            } else {
                if(std::abs(t - mc) > std::abs(t - mq)) {
                    newt = mc;
                } else {
                    newt = mq;
                }
            }
        } else {
            //Case 4: a lower function value, derivatives of the
            //same sign, and the magnitude of the derivative does
            //not decrease. If the minimum is not brackt, the step
            //is either tmin or tmax, else the cubic minimizer is taken.
            bound = 0;
            if(brackt) {
                CUBIC_MINIMIZER(newt, t, ft, dt, y, fy, dy);
            } else if(x < t) {
                newt=tmax;
            } else {
                newt=tmin;
            }
        }

        //Update the interval of uncertainty. This update does not
        //depend on the new step or the case analysis above.
        //
        //	- Case a: if f(x) < f(t),
        //		x <- x, y <- t.
        //	- Case b: if f(t) <= f(x) && f'(t)*f'(x) > 0,
        //		x <- t, y <- y.
        //	- Case c: if f(t) <= f(x) && f'(t)*f'(x) < 0,
        //		x <- t, y <- x.
        if(fx < ft) {
            //Case a
            y=t;
            fy=ft;
            dy=dt;
        } else {
            //Case c
            if(dsign) {
                y=x;
                fy=fx;
                dy=dx;
            }
            //Cases b and c
            x=t;
            fx=ft;
            dx=dt;
        }

        //Clip the new trial value in [tmin, tmax]
        if(tmax < newt)
            newt=tmax;
        if(newt < tmin)
            newt=tmin;

        //Redefine the new trial value if it is close to the upper bound
        //of the interval
        if(brackt && bound) {
            mq=x + 0.66f*(y - x);
            if(x < y) {
                if (mq < newt)
                    newt=mq;
            } else {
                if (newt < mq)
                    newt=mq;
            }
        }

        //Return the new trial value
        t=newt;
        return 0;
    }

#undef USES_MINIMIZER
#undef CUBIC_MINIMIZER
#undef CUBIC_MINIMIZER2
#undef QUARD_MINIMIZER
#undef QUARD_MINIMIZER2
#undef fsigndiff
};

template < typename T, typename KERNEL_TYPE=Kernel<T> >
class LBFGSMinimizer
{
public:
    typedef typename KERNEL_TYPE::Vec Vec;
    struct Correction {
        T _alpha;
        Vec _s;
        Vec _y;
        T _ys;
    };
    //utility
    LBFGSMinimizer() {
        //param settings
        _m=6;
        _epsilon=1e-5f;
        _past=0;
        _delta=1e-5f;
        _maxIterations=0;
        _lineSearch=LINESEARCH_DEFAULT;
        _orthantwiseC=0.0f;
        _orthantwiseStart=0;
        _orthantwiseEnd=0;//-1;

        //no debug
        _probProfile=0;
    }
    sizeType minimize(Vec& x,T& fx,Objective<T,KERNEL_TYPE>& obj,Callback<T,KERNEL_TYPE>& cb=Callback<T,KERNEL_TYPE>()) {
        const sizeType N=obj.inputs();
        ASSERT(x.size() == N)
        //parameter checking
        if(N <= 0)
            return ERR_INVALID_N;
        if(_epsilon < 0.0f)
            return ERR_INVALID_EPSILON;
        if(_past < 0)
            return ERR_INVALID_TESTPERIOD;
        if(_delta < 0.0f)
            return ERR_INVALID_DELTA;
        if(_orthantwiseC < 0.0f)
            return ERR_INVALID_ORTHANTWISE;
        if(_orthantwiseStart < 0 || N < _orthantwiseStart)
            return ERR_INVALID_ORTHANTWISE_START;
        if(_orthantwiseEnd < 0)
            _orthantwiseEnd = N;
        if(N < _orthantwiseEnd)
            return ERR_INVALID_ORTHANTWISE_END;
        if(_orthantwiseC != 0.0f) {
            switch(_lineSearch) {
            case LINESEARCH_BACKTRACKING:
                _searcher.reset(new LineSearcherBackTrackingOWLQN<T>());
                static_cast<LineSearcherBackTrackingOWLQN<T>&>(*_searcher).L1Panelty()=L1Panelty();
                static_cast<LineSearcherBackTrackingOWLQN<T>&>(*_searcher).L1StartIndex()=L1StartIndex();
                static_cast<LineSearcherBackTrackingOWLQN<T>&>(*_searcher).L1EndIndex()=L1EndIndex();
                break;
            default:
                return ERR_INVALID_LINESEARCH;
            }
        } else {
            switch(_lineSearch) {
            case LINESEARCH_MORETHUENTE:
                _searcher.reset(new LineSearcherMoretHuente<T>());
                break;
            case LINESEARCH_BACKTRACKING_ARMIJO:
            case LINESEARCH_BACKTRACKING_WOLFE:
            case LINESEARCH_BACKTRACKING_STRONG_WOLFE:
                _searcher.reset(new LineSearcherBackTracking<T>());
                break;
            default:
                return ERR_INVALID_LINESEARCH;
            }
        }
        _searcher->check(_lineSearch);

        //variables
        fx=0.0f;
        T xnorm=0.0f;
        T gnorm=0.0f;
        T step=0.0f;
        T ys=0.0f;
        T yy=0.0f;
        T beta=0.0f;
        T rate=0.0f;
        sizeType k=0;
        sizeType end=0;
        sizeType bound=0;
        sizeType lsRet=0;
        sizeType cbRet=0;
        Correction* it=NULL;

        //initialization
        _xp.resize(N);	//last x
        _g.resize(N);	//gradient
        _gp.resize(N);	//last gradient
        _d.resize(N);	//search direction
        _w.resize(N);	//variable sign for L1-Opt
        if(_orthantwiseC != 0.0f)
            _pg.resize(N);	//pseudo gradient

        //allocate limited memory updation buffer
        _crr.resize(_m);
        for(sizeType i=0; i<_m; i++) {
            _crr[i]._alpha=0.0f;
            _crr[i]._ys=0.0f;
            _crr[i]._s.resize(N);
            _crr[i]._y.resize(N);
        }

        //Allocate an array for storing previous values of the objective function
        if(0 < _past)
            _pf.resize(_past);

        //compute the function value and derivative
        if(obj(x,fx,_g,0,true) < 0)
            return ERR_USER_ASKED;

        if (0.0f != _orthantwiseC) {
            //Compute the L1 norm of the variable and add it to the object value
            xnorm=L1Norm(x);
            fx+=xnorm*_orthantwiseC;
            L1PseudoGradient(_pg,x,_g);
        }

        //Store the initial value of the objective function
        if(0 < _past)
            _pf[0]=fx;

        //Compute the direction, we assume the initial hessian matrix H_0 as the identity matrix
        if(_orthantwiseC == 0.0f)
            KERNEL_TYPE::ncopy(_g,_d);
        else
            KERNEL_TYPE::ncopy(_pg,_d);

        //Make sure that the initial variables are not a minimizer
        xnorm=KERNEL_TYPE::norm(x);
        if (_orthantwiseC == 0.0f)
            gnorm=KERNEL_TYPE::norm(_g);
        else
            gnorm=KERNEL_TYPE::norm(_pg);
        if (xnorm < 1.0f)
            xnorm=1.0f;
        if (gnorm/xnorm <= _epsilon)
            return ALREADY_MINIMIZED;

        //Compute the initial step:
        //step = 1.0 / sqrt(vecdot(d, d, n))
        step=1.0f/KERNEL_TYPE::norm(_d);

        k=1;
        end=0;
        for(;;) {
            //Store the current position and gradient vectors
            KERNEL_TYPE::copy(x,_xp);
            KERNEL_TYPE::copy(_g,_gp);

            //Search for an optimal step
            if(_probProfile&PROFILE_LINE_SEARCH) {
                obj.profileLineSearch(k,_xp,_d,step);
            }
            if(_orthantwiseC == 0.0f) {
                lsRet=(*_searcher)(x,fx,_g,_d,step,_xp,_gp,_w,obj,_lineSearch);
            } else {
                lsRet=(*_searcher)(x,fx,_g,_d,step,_xp,_gp,_w,obj,_lineSearch);
                L1PseudoGradient(_pg,x,_g);
            }
            if(lsRet < 0) {
                //Revert to the previous point
                KERNEL_TYPE::copy(_xp,x);
                KERNEL_TYPE::copy(_gp,_g);
                return lsRet;
            }

            //Compute x and g norms
            xnorm=KERNEL_TYPE::norm(x);
            if(_orthantwiseC == 0.0f)
                gnorm=KERNEL_TYPE::norm(_g);
            else
                gnorm=KERNEL_TYPE::norm(_pg);

            //Report the progress
            cbRet=cb(x,_g,fx,xnorm,gnorm,step,k,lsRet);
            if(cbRet)
                return cbRet;

            //Convergence test
            //The criterion is given by the following formula:
            //	|g(x)| / \max(1, |x|) < \epsilon
            if (xnorm < 1.0f)
                xnorm=1.0f;
            if (gnorm/xnorm <= _epsilon)
                return SUCCESS;

            //Test for stopping criterion.
            //The criterion is given by the following formula:
            //	(f(past_x) - f(x)) / f(x) < \delta
            if(0 < _past) {
                //We don't test the stopping criterion while k < past
                if(_past <= k) {
                    //Compute the relative improvement from the past
                    rate=(_pf[k%_past]-fx)/fx;
                    //The stopping criterion
                    if (rate < _delta)
                        return STOP;
                }
                //Store the current value of the objective function
                _pf[k%_past]=fx;
            }

            if(_maxIterations != 0 && _maxIterations < k+1)
                return ERR_MAXIMUMITERATION;

            //Update vectors s and y:
            //	s_{k+1} = x_{k+1} - x_{k} = \step * d_{k}.
            //	y_{k+1} = g_{k+1} - g_{k}
            it=&_crr[end];
            KERNEL_TYPE::sub(x,_xp,it->_s);
            KERNEL_TYPE::sub(_g,_gp,it->_y);

            //Compute scalars ys and yy:
            //	ys = y^t \cdot s = 1 / \rho.
            //	yy = y^t \cdot y.
            //Notice that yy is used for scaling the hessian matrix H_0 (Cholesky factor).
            ys=KERNEL_TYPE::dot(it->_y,it->_s);
            yy=KERNEL_TYPE::dot(it->_y,it->_y);
            it->_ys=ys;

            //Recursive formula to compute dir = -(H \cdot g).
            //	This is described in page 779 of:
            //	Jorge Nocedal.
            //	Updating Quasi-Newton Matrices with Limited Storage.
            //	Mathematics of Computation, Vol. 35, No. 151,
            //	pp. 773--782, 1980.
            bound=(_m <= k) ? _m : k;
            ++k;
            end=(end+1)%_m;

            //Compute the steepest direction
            if(_orthantwiseC == 0.0f)
                //Compute the negative of gradients
                KERNEL_TYPE::ncopy(_g,_d);
            else
                KERNEL_TYPE::ncopy(_pg,_d);

            sizeType j=end;
            for (sizeType i=0; i < bound; ++i) {
                j=(j+_m-1)%_m;
                it=&_crr[j];
                it->_alpha=KERNEL_TYPE::dot(it->_s,_d);
                it->_alpha/=it->_ys;
                KERNEL_TYPE::addScaled(-it->_alpha,it->_y,_d);
            }
            KERNEL_TYPE::scale(ys/yy,_d);
            for (sizeType i=0; i < bound; ++i) {
                it=&_crr[j];
                beta=KERNEL_TYPE::dot(it->_y,_d);
                beta/=it->_ys;
                KERNEL_TYPE::addScaled(it->_alpha-beta,it->_s,_d);
                j=(j+1)%_m;
            }

            //Constrain the search direction for orthant-wise updates
            if(_orthantwiseC != 0.0f) {
                for(sizeType i=_orthantwiseStart; i < _orthantwiseEnd; ++i) {
                    if(_d[i]*_pg[i] >= 0.0f)
                        _d[i]=0.0f;
                }
            }

            //Now the search direction d is ready. We try step = 1 first
            step=1.0f;
        }

        return SUCCESS;
    }
    //param
    const sizeType& nrCorrect() const {
        return _m;
    }
    sizeType& nrCorrect() {
        return _m;
    }
    const T& epsilon() const {
        return _epsilon;
    }
    T& epsilon() {
        return _epsilon;
    }
    const sizeType& historyLengthOfConvergenceTest() const {
        return _past;
    }
    sizeType& historyLengthOfConvergenceTest() {
        return _past;
    }
    const T& delta() const {
        return _delta;
    }
    T& delta() {
        return _delta;
    }
    const sizeType& maxIterations() const {
        return _maxIterations;
    }
    sizeType& maxIterations() {
        return _maxIterations;
    }
    const LINE_SEARCHER& lineSearch() const {
        return _lineSearch;
    }
    LINE_SEARCHER& lineSearch() {
        return _lineSearch;
    }
    const T& L1Panelty() const {
        return _orthantwiseC;
    }
    T& L1Panelty() {
        return _orthantwiseC;
    }
    const sizeType& L1StartIndex() const {
        return _orthantwiseStart;
    }
    sizeType& L1StartIndex() {
        return _orthantwiseStart;
    }
    const sizeType& L1EndIndex() const {
        return _orthantwiseEnd;
    }
    sizeType& L1EndIndex() {
        return _orthantwiseEnd;
    }
    const sizeType& probProfile() const {
        return _probProfile;
    }
    sizeType& probProfile() {
        return _probProfile;
    }
protected:
    //L1-norm minimizer helper
    T L1Norm(const Vec& x) const {
        T norm=0.0f;
        for(sizeType i=_orthantwiseStart; i<_orthantwiseEnd; ++i)
            norm+=std::abs(x[i]);
        return norm;
    }
    void L1PseudoGradient(Vec& pg,const Vec& x,const Vec& g) const {
        const sizeType N=x.size();
        for(sizeType i=0; i<_orthantwiseStart; ++i)
            pg[i]=g[i];
        for(sizeType i=_orthantwiseStart; i<_orthantwiseEnd; ++i) {
            if(x[i] < 0.0f)
                pg[i]=g[i]-_orthantwiseC;
            else if(0.0f < x[i])
                pg[i]=g[i]+_orthantwiseC;
            else {
                if(g[i] < -_orthantwiseC)
                    pg[i]=g[i]+_orthantwiseC;
                else if(_orthantwiseC < g[i])
                    pg[i]=g[i]-_orthantwiseC;
                else
                    pg[i]=0.0f;
            }
        }
        for(sizeType i=_orthantwiseEnd; i<N; ++i)
            pg[i]=g[i];
    }
    //params
    sizeType _m;
    T _epsilon;
    sizeType _past;
    T _delta;
    sizeType _maxIterations;
    LINE_SEARCHER _lineSearch;
    T _orthantwiseC;
    sizeType _orthantwiseStart;
    sizeType _orthantwiseEnd;
    //data
    boost::shared_ptr<LineSearcher<T> > _searcher;
    Vec _xp,_g,_pg,_gp,_d,_w,_pf;
    std::vector<Correction> _crr;
    //debugger visualization
    sizeType _probProfile;
};

PRJ_END

#endif
