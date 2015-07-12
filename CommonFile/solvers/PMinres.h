#ifndef PMINRES_H
#define PMINRES_H

#include "LinearSolver.h"
#include <iostream>
#include <iomanip>

PRJ_BEGIN

template <typename T,typename KERNEL_TYPE=Kernel<T> >
class KKTKrylovMatrix : public KrylovMatrix<T,KERNEL_TYPE>
{
public:
    typedef typename KrylovMatrix<T,KERNEL_TYPE>::Vec Vec;
	KKTKrylovMatrix(const FixedSparseMatrix<T,KERNEL_TYPE>& A,
					const FixedSparseMatrix<T,KERNEL_TYPE>& C):_A(A),_C(C){}
	virtual void multiply(const Vec& b,Vec& out) const
	{
		Eigen::Block<Vec> outVec=out.block(_A.rows(),0,_C.rows(),1);
		KERNEL_TYPE::zero(out);
		_A.multiplyAdd(b,out);
		_C.multiplyAdd(b,outVec);
		_C.multiplyTransposeAdd(b.block(_A.rows(),0,_C.rows(),1),out);
	}
	virtual sizeType n() const{return _A.rows()+_C.rows();}
	//data
	const FixedSparseMatrix<T,KERNEL_TYPE>& _A;
	const FixedSparseMatrix<T,KERNEL_TYPE>& _C;
};
template <typename T,typename KERNEL_TYPE=Kernel<T> >
class EigenKKTKrylovMatrix : public KrylovMatrix<T,KERNEL_TYPE>
{
public:
    typedef typename KrylovMatrix<T,KERNEL_TYPE>::Vec Vec;
	EigenKKTKrylovMatrix(const Eigen::SparseMatrix<T,0,sizeType>& A,const Eigen::SparseMatrix<T,0,sizeType>& C):_A(A),_C(C){}
	virtual void multiply(const Vec& b,Vec& out) const
	{
		out.block(0,0,_A.rows(),1)=_A*b.block(0,0,_A.rows(),1);
		out.block(0,0,_A.rows(),1)+=_C.transpose()*b.block(_A.rows(),0,_C.rows(),1);
		out.block(_A.rows(),0,_C.rows(),1)=_C*b.block(0,0,_A.rows(),1);
	}
	virtual sizeType n() const{return _A.rows()+_C.rows();}
	//data
	const Eigen::SparseMatrix<T,0,sizeType>& _A;
	const Eigen::SparseMatrix<T,0,sizeType>& _C;
};

//============================================================================
// PMINRES: Just a port of Choi's code (stanford SOL lib)
template <typename T,typename KERNEL_TYPE=Kernel<T>,typename PRECONDITIONER_TYPE=NoPreconSolver<T> >
struct PMINRESSolver : public Solver<T,KERNEL_TYPE> {
public:
    typedef typename KERNEL_TYPE::Vec Vec;
    using Solver<T,KERNEL_TYPE>::_cb;
    using Solver<T,KERNEL_TYPE>::_residualOut;
    using Solver<T,KERNEL_TYPE>::_iterationsOut;
    PMINRESSolver(void) {
        setSolverParameters(1e-5f,1000);
        _debug=false;
        _shift=0.0f;
    }
    void setDebug(bool debug){
        _debug=debug;
    }
    void setMaxxnorm(T maxxnorm){
        _maxxnorm=maxxnorm;
    }
    void setACondlim(T Acondlim){
        _Acondlim=Acondlim;
    }
    virtual void setSolverParameters(T toleranceFactor,sizeType maxIterations) {
        _toleranceFactor=toleranceFactor;
        if(_toleranceFactor<1e-30f)
            _toleranceFactor=1e-30f;
        _maxxnorm=1E14f;
        _Acondlim=1E14f;
        _maxIterations=maxIterations;
    }
    virtual void resize(const sizeType& n) {
        if(_r1.size()!=n) {
            _r1.resize(n);
            _r2.resize(n);
            _r3.resize(n);
            _d.resize(n);
            _dl.resize(n);
            _dl2.resize(n);
            _v.resize(n);
        }
    }
    virtual void setMatrix(const FixedSparseMatrix<T,KERNEL_TYPE>& matrix,bool syncPrecon) {
        resize(matrix.rows());
        if(syncPrecon)_pre.setMatrix(matrix,true);
        _fixedMatrix.reset(new DefaultKrylovMatrix<T,KERNEL_TYPE>(matrix));
    }
    virtual void setKrylovMatrix(boost::shared_ptr<KrylovMatrix<T,KERNEL_TYPE> > matrix) {
        resize(matrix->n());
        _fixedMatrix=matrix;
    }
    virtual typename Solver<T,KERNEL_TYPE>::SOLVER_RESULT solve(const Vec& rhs,Vec& x) {
        // Initalize
        sizeType flag=0;
        bool done=false;
        
        std::vector<std::string> msg(13);
        msg[0]=" beta2 = 0.  If M = I, b and x are eigenvectors.              ";
        msg[1]=" beta1 = 0.  The exact solution is  x = 0.                    ";   //  0
        msg[2]=" A solution to (poss. singular) Ax = b found, given rtol.     ";   //  1
        msg[3]=" A least-squares solution was found, given rtol.              ";   //  2
        msg[4]=" A solution to (poss. singular) Ax = b found, given eps.      ";   //  3
        msg[5]=" A least-squares solution was found, given eps.               ";   //  4
        msg[6]=" x has converged to an eigenvector.                           ";   //  5
        msg[7]=" xnorm has exceeded maxxnorm.                                 ";   //  6
        msg[8]=" Acond has exceeded Acondlim.                                 ";   //  7
        msg[9]=" The iteration limit was reached.                             ";   //  8
       msg[10]=" A least-squares solution for singular LS problem, given eps. ";   //  9
       msg[11]=" A least-squares solution for singular LS problem, given rtol.";   //  10
       msg[12]=" A null vector obtained, given rtol.                          ";   //  11

        T beta1=sqrt(KERNEL_TYPE::dot(rhs,rhs));
        if(_debug){
            INFOV("shift=%f, tol=%f, maxIter=%d",_shift,_toleranceFactor,(int)_maxIterations)
            std::cout<<std::setfill('-')<<std::setw(80)<<"-"<< "\n";
            std::cout<<"|            Adapted from minresSOL69.m, Stanford University, 06 Jul 2009       |\n";
            std::cout<<"|                Solution of symmetric Ax=b or (A-shift*I)x=b                   |\n";
            std::cout<<std::setfill('-')<<std::setw(80)<<"-"<<"\n";
            std::cout<<std::setfill(' ');
        }
        if(_cb){
            _cb->reset();
        }

        //------------------------------------------------------------------
        // Set up p and v for the first Lanczos vector v1.
        // p=beta1 P' v1,  where  P=C**(-1).
        // v is really P' v1.
        //------------------------------------------------------------------
        KERNEL_TYPE::copy(rhs,_r3);
        KERNEL_TYPE::copy(rhs,_r2);
        KERNEL_TYPE::copy(rhs,_r1);

        _pre.solve(rhs,_r3);
        beta1=KERNEL_TYPE::dot(_r3,rhs);

        //  Test for an indefinite preconditioner.
        //  If b=0 exactly, stop with x=0.
        if(beta1 < 0.0f){
            WARNING("The preconditioner 'M' appears to be indefinite");
        }else{
            beta1=sqrt(beta1);
        }

        if(beta1 == 0.0f){ // b=0, so x=0
            flag=0;
            done=true;
        }

        //------------------------------------------------------------------
        // Initialize other quantities.
        // Note that Anorm has been initialized by IsOpSym6.
        //------------------------------------------------------------------
        sizeType& iter=_iterationsOut=0;
        T beta=0.0f;
        T phi=beta1;
        T Acond=1.0f;
        T betan=beta1;
        T bnorm=beta1;
        T cs=-1.0f;
        T dbarn=0.0f;
        T eplnn=0.0f;
        T rnorm=betan;
        T sn=0.0f;
        T Tnorm=0.0f;
        T rnorml=beta1;
        T xnorm=0.0f;
        T Dnorm=0.0f;
        //T xnorm=0.0f;
        T xrArnormsl=0.0f;
        T gama=0.0f;
        sizeType lastiter=0;
        T pnorm=0.0f;
        T gamal=0.0f;
        T Axnorm=0.0f;
        T relrnorm=rnorm/beta1;
        T Anorm=0.0f;

        T eps=1E-14f;//ScalarUtil<T>::scalar_eps;
        T betal,pnorml,xnorml,x1last,rootl,relrnorml,relArnorml,Arnorml,Anorml,Acondl;
        T Arnorm,relArnorm,d_norm,alfa,dbar,epln,dlta,gbar,gamal2,tau,t2,t1,epsx,epsr;

        _fixedMatrix->multiply(rhs,_dl);  // dl temporarily used for testing
        Arnorm=sqrt(KERNEL_TYPE::dot(_dl,_dl));

        KERNEL_TYPE::copy(rhs,x);
        KERNEL_TYPE::zero(x);
        KERNEL_TYPE::zero(_d);
        KERNEL_TYPE::zero(_dl);

        if(_debug){
            std::cout<<std::setw(6)<<"Iter"
                     << std::setw(14) << "x(1)"
                     << std::setw(14) << "xnorm"
                     << std::setw(14) << "rnorm"
                     << std::setw(14) << "Arnorm"
                     << std::setw(14) << "Compatible"
                     << std::setw(14) << "LS"
                     << std::setw(14) << "norm(A)"
                     << std::setw(14) << "cond(A)"
                     << std::setw(14) << "beta"<<"\n";
        }

        //---------------------------------------------------------------------
        // Main iteration loop.
        //--------------------------------------------------------------------
        while(!done){
            iter=iter+1;
            //-----------------------------------------------------------------
            // Obtain quantities for the next Lanczos vector vk+1, k=1, 2,...
            // The general iteration is similar to the case k=1 with v0=0:
            //
            //   p1=Operator * v1  -  beta1 * v0,
            //   alpha1=v1'p1,
            //   q2=p2  -  alpha1 * v1,
            //   beta2^2=q2'q2,
            //   v2=(1/beta2) q2.
            //
            // Again, p=betak P vk,  where  P=C**(-1).
            // .... more description needed.
            //-----------------------------------------------------------------
            betal=beta;
            beta=betan;
            KERNEL_TYPE::copy(_r3,_v);
            KERNEL_TYPE::scale(1.0f/beta,_v);                   //v=r3 * (1/beta);
            _fixedMatrix->multiply(_v,_r3);                     // p_k
            if(_shift != 0.0f){
                KERNEL_TYPE::addScaled(-_shift,_v,_r3);         //r3=r3 - shift*v;
            }

            if(iter > 1){
                KERNEL_TYPE::addScaled(-(beta/betal),_r1,_r3);  //r3=r3 - (beta/betal) * r1;
            }
            alfa=KERNEL_TYPE::dot(_r3,_v);                      //alfa=real( r3'*v );                           // alpha_k
            KERNEL_TYPE::addScaled(-(alfa/beta),_r2,_r3);       //r3=r3 - (alfa/beta ) * r2;                    // v_{k+1}
            KERNEL_TYPE::copy(_r2,_r1);                         //r1=r2;                                        // v_{k-1}
            KERNEL_TYPE::copy(_r3,_r2);                         //r2=r3;

            _pre.solve(_r2,_r3);
            betan=KERNEL_TYPE::dot(_r2,_r3);                    //betan=r2' * r3;
            if(betan < 0.0f){
                WARNING("The preconditioner 'M' appears to be indefinite");
            }
            if(betan > 0.0f){
                betan=sqrt(betan);
            }

            pnorml=pnorm;
            if(iter == 1){
                pnorm=sqrt(alfa*alfa+betan*betan);
            } else {
                pnorm=sqrt(beta*beta+alfa*alfa+betan*betan);
            }

            //-----------------------------------------------------------------
            // Apply previous rotation Qk-1 to get
            //   [dlta_k epln_{k+1}]=[cs  sn][dbar_k    0      ]
            //   [gbar_k  dbar_{k+1} ]   [sn -cs][alfa_k beta_{k+1}].
            //-----------------------------------------------------------------
            dbar=dbarn;
            epln=eplnn;

            dlta=cs*dbar+sn*alfa;               // dlta1=0         dltak
            gbar=sn*dbar-cs*alfa;               // gbar 1=alfa1     gbar k

            eplnn=sn*betan;                     // epln2=0        eplnk+1
            dbarn=-cs*betan;                    // dbarn 2=beta2   dbarn k+1

            // Compute the current plane rotation Qk
            gamal2=gamal;
            gamal=gama;
            symGivens2(gbar,betan,  cs,sn,gama);
            tau=cs*phi;                         // phik
            phi=sn*phi;                         // phibar_{k+1}
            Axnorm=sqrt(Axnorm*Axnorm+tau*tau); // norm([tau_1,..., tau_k])

            // Update d
            KERNEL_TYPE::copy(_dl,_dl2);        //dl2=dl;
            KERNEL_TYPE::copy(_d,_dl);          //dl=d;
            if(gama !=0.0f) {
                KERNEL_TYPE::copy(_v,_d);
                KERNEL_TYPE::addScaled(-epln,_dl2,_d);
                KERNEL_TYPE::addScaled(-dlta,_dl,_d);
                KERNEL_TYPE::scale(1.0f/gama,_d);       //d=(v - epln*dl2 - dlta*dl) * (1/gama);
                d_norm=sqrt(KERNEL_TYPE::dot(_d,_d));   //d_norm=norm(d);
            } else {
                d_norm=1.0f/eps;                        //d_norm=Inf;
            }

            // Update x except if it will become too big
            x1last=x(0);
            xnorml=xnorm;
            KERNEL_TYPE::copy(x,_dl2);                  //dl2=x;   // temporarily using dl2 to store x_{k-1}
            KERNEL_TYPE::addScaled(tau,_d,x);           //x=x+tau*d;
            xnorm=sqrt(KERNEL_TYPE::dot(x,x));          //xnorm=norm(x);
            if(xnorm >=_maxxnorm){                      //max(d_norm,xnorm) >=_maxxnorm
                KERNEL_TYPE::copy(_dl2,x);              //x=dl2;
                flag=6;
                lastiter=1;
            }

            // Estimate various norms
            //xnorml=xnorm;
            rnorml=rnorm;                               // ||r_{k-1}||
            Anorml=Anorm;
            Acondl=Acond;
            relrnorml=relrnorm;

            if(flag != 6){
                Dnorm=sqrt(Dnorm*Dnorm+d_norm*d_norm);
                xnorm=sqrt(KERNEL_TYPE::dot(x,x));      // rather expensive

                rnorm=phi;                              // ||r_k||
                relrnorm=rnorm/(Anorm*xnorm+bnorm);     // ||r||/(||A||||x||+||b||)

                if(iter==1) {
                    Tnorm=sqrt(alfa*alfa+betan*betan);
                } else {
                    Tnorm=sqrt(Tnorm*Tnorm+beta*beta+alfa*alfa+betan*betan);
                    // Frobenius norm of T_{k,k+1}
                }   // used as estimate of Anorm.
                //Anorm=Tnorm;
                Anorm=std::max(Anorm,pnorm);
                Acond=Anorm*Dnorm;
            }

            rootl=sqrt(gbar*gbar+dbarn*dbarn);
            Arnorml=rnorml*rootl;                       // ||A r_{k-1} ||
            relArnorml=rootl/Anorm;                     // ||Ar|| / (||A|| ||r||)
            //relArnorml=Arnorml/Anorm;                 // ||Ar|| / ||A||

            //---------------------------------------------------------------
            // See if any of the stopping criteria are satisfied.
            // In rare cases, flag is already -1 from above (Abar=const*I).
            //---------------------------------------------------------------
            // epsa=Anorm * eps;                        // Used below to test betan and ||A r||
            epsx=Anorm*xnorm*eps;
            epsr=Anorm*xnorm*_toleranceFactor;

            // Test for singular Hk (hence singular A)
            // or x is already an LS solution (so again A must be singular).
            if(flag==0 || flag==6) {
                t1=1.0f+relrnorm;                                   // These tests work if _toleranceFactor < eps
                t2=1.0f+relArnorml;
                if(_toleranceFactor >= eps) {
                    if(t1 <= 1.0f)flag=3;                            // Accurate Ax=b solution
                    else if(t2 <= 1.0f)flag=4;                       // Accurate LS solution
                    else if(relrnorm <= _toleranceFactor)flag=1;     // Good enough Ax=b solution
                    else if(relArnorml <= _toleranceFactor)flag=10;  // Good enough LS solution
                    else if(epsx >= beta1)flag=5;                    // x is an eigenvector
                    else if(xnorm >= _maxxnorm)flag=6;               // xnorm exceeded its limit
                    else if(Acond >= _Acondlim)flag=7;               // Huge Acond (but maybe ok?)
                    else if(iter >= _maxIterations)flag=8;           // Too many itns
                } else {
                    if(relrnorm <= _toleranceFactor)flag=1;          // Good enough Ax=b solution
                    else if(relArnorml <= _toleranceFactor)flag=10;  // Good enough LS solution
                    else if(epsx >= beta1)flag=5;                    // x is an eigenvector
                    else if(xnorm >= _maxxnorm)flag=6;               // xnorm exceeded its limit
                    else if(Acond >= _Acondlim)flag=7;               // Huge Acond (but maybe ok?)
                    else if(iter >= _maxIterations)flag=8;           // Too many itns
                }
            }
            if(flag != 0){                                                  // flag may have been changed above
                done=1;
                if((flag==1) || (flag==3) || (flag==8) || (flag==11)) {     // Linear systems
                    _fixedMatrix->multiply(x,_r1);
                    KERNEL_TYPE::scale(-1.0f,_r1);
                    KERNEL_TYPE::add(_r1,rhs,_r1);
                    KERNEL_TYPE::addScaled(_shift,x,_r1);                   //r1=b-minresxxxA( A, x , varargin )+_shift*x;      //use r1 to temporarily store residual
                    
                    _fixedMatrix->multiply(_r1,_r2);
                    KERNEL_TYPE::addScaled(-_shift,_r1,_r2);
                    Arnorm=sqrt(KERNEL_TYPE::dot(_r2,_r2));                 //Arnorm=norm(minresxxxA(A,r1,varargin )-_shift*r1); // avoid computing 1 more Lanczos iteration to get root
                    
                    relArnorm=Arnorm/(Anorm*rnorm);                         // ||Ar|| / ||A||
                }
                if((flag==4) || (flag==6) ||  (flag==10) || (flag==7)) {
                    lastiter=1;
                }
                if(lastiter==1) {
                    iter=iter-1;
                    rnorm=rnorml;
                    relrnorm=relrnorml;
                    Arnorm=Arnorml;
                    relArnorm=relArnorml;
                    pnorm=pnorml;
                    gama=gamal;
                    gamal=gamal2;
                }
            }

            // See if it is time to print something.
            if(_debug) {
                std::cout<< std::setw(6) << iter-1
                         << std::setw(14) << x1last
                         << std::setw(14) << xnorml
                         << std::setw(14) << rnorml
                         << std::setw(14) << Arnorml
                         << std::setw(14) << relrnorml
                         << std::setw(14) << relArnorml
                         << std::setw(14) << Anorml
                         << std::setw(14) << Acondl
                         << std::setw(14) << betal << std::endl;
            }
            _residualOut=relrnorml;
            if(_cb){
                if((*_cb)(x,relrnorml,iter))
                    return Solver<T,KERNEL_TYPE>::USER_REQUEST;
            }
        }

        if(_debug) {
            std::cout<<std::setfill('-')<<std::setw(80)<<"-"<< "\n";
            std::cout<<"|                                  Finish                                       |\n";
            std::cout<<std::setfill('-')<<std::setw(80)<<"-"<< "\n";
            std::cout<<std::setfill(' ');

            std::cout<<std::setw(6)<<"Iter"
                     << std::setw(14) << "x(1)"
                     << std::setw(14) << "xnorm"
                     << std::setw(14) << "rnorm"
                     << std::setw(14) << "Arnorm"
                     << std::setw(14) << "Compatible"
                     << std::setw(14) << "LS"
                     << std::setw(14) << "norm(A)"
                     << std::setw(14) << "cond(A)"
                     << std::setw(14) << "beta"<<"\n";
            std::cout<< std::setw(6) << iter-1
                     << std::setw(14) << x1last
                     << std::setw(14) << xnorml
                     << std::setw(14) << rnorml
                     << std::setw(14) << Arnorml
                     << std::setw(14) << relrnorml
                     << std::setw(14) << relArnorml
                     << std::setw(14) << Anorml
                     << std::setw(14) << Acondl
                     << std::setw(14) << betal << std::endl;
			std::cout << msg[flag+1] << std::endl;
        }

        return (flag == 6 || flag == 7 || flag == 8) ? Solver<T,KERNEL_TYPE>::NOT_CONVERGENT : Solver<T,KERNEL_TYPE>::SUCCESSFUL;
    }
    Solver<T,KERNEL_TYPE>* getPre() {
        return &_pre;
    }
    const Solver<T,KERNEL_TYPE>* getPre() const {
        return &_pre;
    }
    static void symGivens2(T a,T b, T& c,T& s,T& d) {
        if(b==0.0f) {
            if(a==0.0f) {
                c=1.0f;
            } else {
                c=(T)sgn(a);  // NOTE: sign(0)=0 in MATLAB
            }
            s=0.0f;
            d=std::abs(a);
        } else if(a==0.0f) {
            c=0.0f;
            s=(T)sgn(b);
            d=std::abs(b);
        } else if(std::abs(b) > std::abs(a)) {
            T t=a/b;
            s=(T)sgn(b)/sqrt(1.0f+t*t);
            c=s*t;
            d=b/s;          // computationally better than d=a / c since |c| <=|s|
        } else {
            T t=b/a;
            c=(T)sgn(a)/sqrt(1.0f+t*t);
            s=c*t;
            d=a/c;          // computationally better than d=b / s since |s| <=|c|
        }
        // end of minresQLP
    }
protected:
    // internal structures
    PRECONDITIONER_TYPE _pre;
    Vec _r1,_r2,_r3,_d,_dl,_dl2,_v; // temporary vectors for PCG
    boost::shared_ptr<KrylovMatrix<T,KERNEL_TYPE> > _fixedMatrix;
    // parameters
    T _toleranceFactor;
    T _shift;
    T _maxxnorm,_Acondlim;
    sizeType _maxIterations;
    bool _debug;
};

PRJ_END

#endif
