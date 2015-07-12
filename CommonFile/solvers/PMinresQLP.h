#ifndef PMINRES_QLP_H
#define PMINRES_QLP_H

#include "PMinres.h"

PRJ_BEGIN

template <typename T,typename KERNEL_TYPE=Kernel<T>,typename PRECONDITIONER_TYPE=NoPreconSolver<T> >
struct PMINRESSolverQLP : public Solver<T,KERNEL_TYPE> {
public:
    typedef typename KERNEL_TYPE::Vec Vec;
    using Solver<T,KERNEL_TYPE>::_cb;
    using Solver<T,KERNEL_TYPE>::_residualOut;
    using Solver<T,KERNEL_TYPE>::_iterationsOut;
    PMINRESSolverQLP(void) {
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
            _xl2.resize(n);
            _r1.resize(n);
            _r2.resize(n);
            _r3.resize(n);
            _w.resize(n);
            _wl.resize(n);
            _wl2.resize(n);
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
        sizeType flag=0;
        bool done=false;

        std::vector<std::string> msg(12);
        msg[0]=" beta2 = 0.  If M = I, b and x are eigenvectors.             ";   // -1
        msg[1]=" beta1 = 0.  The exact solution is  x = 0.                   ";   //  0
        msg[2]=" A solution to (poss. singular) Ax = b found, given rtol.    ";   //  1
        msg[3]=" Pseudoinverse solution for singular LS problem, given rtol. ";   //  2
        msg[4]=" A solution to (poss. singular) Ax = b found, given eps.     ";   //  3
        msg[5]=" Pseudoinverse solution for singular LS problem, given eps.  ";   //  4
        msg[6]=" x has converged to an eigenvector.                          ";   //  5
        msg[7]=" xnorm has exceeded maxxnorm.                                ";   //  6
        msg[8]=" Acond has exceeded Acondlim.                                ";   //  7
        msg[9]=" The iteration limit was reached.                            ";   //  8
       msg[10]=" It is a least squares problem but no converged solution yet.";   //  9
       msg[11]=" A null vector obtained, given rtol.                         ";  //  10

        T beta1=sqrt(KERNEL_TYPE::dot(rhs,rhs));
        if(_debug) {
            INFOV("shift=%f, tol=%f, maxIter=%d",_shift,_toleranceFactor,_maxIterations)
            std::cout<<std::setfill('-')<<std::setw(80)<<"-"<< "\n";
            std::cout<<"|            Adapted from minresQLP49.m, Stanford University, 06 Jul 2009       |\n";
            std::cout<<"|            Minimum Length Solution of symmetric Ax=b or (A-shift*I)x=b        |\n";
            std::cout<<std::setfill('-')<<std::setw(80)<<"-"<<"\n";
            std::cout<<std::setfill(' ');
        }
        if(_cb) {
            _cb->reset();
        }

        //------------------------------------------------------------------
        // Set up p and v for the first Lanczos vector v1.
        // p  =  beta1 P' v1,  where  P = inv(C).
        // v is really P' v1.
        //------------------------------------------------------------------
        KERNEL_TYPE::copy(rhs,_r3);//r3    = b;
        KERNEL_TYPE::copy(rhs,_r2);//r2    = b;
        _pre.solve(rhs,_r3);
        beta1=KERNEL_TYPE::dot(_r3,rhs);
        //  Test for an indefinite preconditioner.
        //  If b = 0 exactly, stop with x = 0.
        if(beta1 >= 0.0f) {
            beta1=sqrt(beta1);
        } else {
            WARNINGV("The preconditioner 'M' appears to be indefinite or singular.");
        }
        if(beta1 == 0.0f) {
            // b = 0, so x = 0
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
        T betan=beta1;
        T gmin=0.0f;
        T cs=-1.0f;
        T sn=0.0f;
        T cr1=-1.0f;
        T sr1=0.0f;
        T cr2=-1.0f;
        T sr2=0.0f;
        T dltan=0.0f;
        T gama=0.0f;
        T gamal=0.0f;
        T eta=0.0f;
        T etal=0.0f;
        T etal2=0.0f;
        T vepln=0.0f;
        T veplnl=0.0f;
        T veplnl2=0.0f;
        T taul=0.0f;
        T tau=0.0f;
        T Acond=1.0f;
        T rnorm=betan;
        T pnorm=0.0f;
        T xnorm=0.0f;
        T Anorm=0.0f;
        T xl2norm=0.0f;
        T ul3=0.0f;
        T ul2=0.0f;
        T ul=0.0f;
        sizeType lastiter=0;
        //bestiter=0;
        //prnt2=true;
        T Axnorm=0.0f;
        //lines=1;
        //headlines=30 * lines;
        T L_ek_norm=0.0f;
        T u=0.0f;
        T gamal2=0.0f;

        T relres=rnorm/(Anorm*xnorm+beta1);
        //rnormbest=rnorm;
        //BIG=1000;
        //xrArnormsbest=BIG*Anorm*rnorm^2;

        T xnorml,xnorm_tmp,Tnorm,rnorml,pnorml,Anorml,Arnorml,L_el2_norm,L_ek_norml,Acondl;//,AnormF,relrnorml
        T relresl,relAres,relAresl,x1last,ul4,taul2,t1,t2,rootl,gbar,gamal3,gamal_tmp,gama_tmp;
        T epsx,epsr,eplnn,dlta_tmp,dlta,dbar,betal,alfa;
        T eps=ScalarUtil<T>::scalar_eps;

        _fixedMatrix->multiply(rhs,_w);

        T Arnorm=sqrt(KERNEL_TYPE::dot(_w,_w));

        KERNEL_TYPE::zero(_xl2);
        KERNEL_TYPE::copy(rhs,x);
        KERNEL_TYPE::zero(x);
        KERNEL_TYPE::zero(_w); // reset w to 0
        KERNEL_TYPE::zero(_wl);

        if(_debug) {
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

        while(!done && iter < _maxIterations) {
            iter=iter+1;

            //-----------------------------------------------------------------
            // Obtain quantities for the next Lanczos vector vk+1, k = 1, 2,...
            // The general iteration is similar to the case k = 1 with v0 = 0:
            //
            //   p1      = Operator * v1  -  beta1 * v0,
            //   alpha1  = v1'p1,
            //   q2      = p2  -  alpha1 * v1,
            //   beta2^2 = q2'q2,
            //   v2      = (1/beta2) q2.
            //
            // Again, p = betak P vk,  where  P = C**(-1).
            // .... more description needed.
            //-----------------------------------------------------------------

            // Lanczos
            betal=beta;
            beta=betan;
            KERNEL_TYPE::copy(_r3,_v);
            KERNEL_TYPE::scale(1.0f/beta,_v);                       //v=r3*(1/beta);// If precon, v=s_k / beta k.
            _fixedMatrix->multiply(_v,_r3);
            if(_shift != 0.0f) {
                KERNEL_TYPE::addScaled(-_shift,_v,_r3);             //r3=r3 - shift*v;// r3=p_k (will be changed below).//If precon, r3=p_k / beta_k
            }

            if(iter > 1) {
                KERNEL_TYPE::addScaled(-(beta/betal),_r1,_r3);      //r3=r3 - (beta/betal) * r1;// r1=beta_{k-1} v_{k-1}. If precon, r1=z_k-1.
            }
            alfa=KERNEL_TYPE::dot(_r3,_v);                          //alfa=real( r3' * v );// alpha_k
            KERNEL_TYPE::addScaled(-(alfa/beta),_r2,_r3);           //r3=r3 - ( alfa/beta ) * r2;// r2=beta_k v_k, r3=beta_{k+1} v_{k+1}.// If precon, r2=z_k, r3=z_k+1.
            KERNEL_TYPE::copy(_r2,_r1);        //r1=r2;// If precon, r1=z_k
            KERNEL_TYPE::copy(_r3,_r2);        //r2=r3;// If precon, r2=z_k+1
            {
                _pre.solve(_r2,_r3);
                betan=KERNEL_TYPE::dot(_r2,_r3);    //betan=r2' * r3;
                if(betan >= 0.0f) {
                    betan=sqrt(betan);
                } else {
                    WARNINGV("The preconditioner 'M' appears to be indefinite or singular.");
                }
            }
            pnorml=pnorm;
            if(iter==1) {
                pnorm=sqrt(alfa*alfa+betan*betan);
            } else {
                pnorm=sqrt(beta*beta+alfa*alfa+betan*betan);
            }

            //-----------------------------------------------------------------
            // Apply previous left rotation Qk-1 to get
            //   [dlta_k epln_{k+1} ] = [cs  sn][dlta_k    0      ]
            //   [gbar_k  dlta_{k+1}]   [sn -cs][alfa_k beta_{k+1}].
            //-----------------------------------------------------------------
            dbar=dltan;
            dlta=cs*dbar+sn*alfa;       // dlta1=0         dlta k
            gbar=sn*dbar-cs*alfa;       // gbar 1=alfa1     gbar k
            eplnn=sn*betan;             // epln2=0         epln k+1
            dltan=-cs*betan;            // dltan 2=beta2   dltan k+1

            // Compute the current left plane rotation Qk
            gamal3=gamal2;
            gamal2=gamal;
            gamal=gama;
            PMINRESSolver<T>::symGivens2(gbar,betan, cs,sn,gama);
            gama_tmp=gama;
            taul2=taul;
            taul=tau;
            tau=cs*phi;                             // phik
            phi=sn*phi;                             // phibar_{k+1}
            Axnorm=sqrt(Axnorm*Axnorm+tau*tau);     // norm([tau_1,..., tau_k])

            //  Apply the previous right plane rotation P{k-2,k}
            if(iter > 2) {
                veplnl2=veplnl;
                etal2=etal;
                etal=eta;
                dlta_tmp=sr2*vepln-cr2*dlta;        // vepln=vepln from last iter
                veplnl=cr2*vepln+sr2*dlta;
                dlta=dlta_tmp;
                eta=sr2*gama;
                gama=-cr2*gama;
            }

            // Compute the current right plane rotation P{k-1,k}, P_12, P_23,...
            if(iter > 1) {
                PMINRESSolver<T>::symGivens2(gamal,dlta, cr1,sr1,gamal);
                vepln=sr1*gama;
                gama=-cr1*gama;
            }

            // Update xnorm
            xnorml=xnorm;
            ul4=ul3;
            ul3=ul2;
            if(iter > 2) {
                ul2=(taul2-etal2*ul4-veplnl2*ul3)/gamal2;
            }
            if(iter > 1) {
                ul=(taul-etal*ul3-veplnl*ul2)/gamal;
            }
            xnorm_tmp=sqrt(xl2norm*xl2norm+ul2*ul2+ul*ul);
            if((std::abs(gama) > eps) && (xnorm_tmp < _maxxnorm)) {
                u=(tau-eta*ul2-vepln*ul)/gama;
                if(sqrt(xnorm_tmp*xnorm_tmp+u*u) > _maxxnorm) {
                    u=0;
                    flag=6;
                }
            } else {
                u=0;
                flag=9;
            }
            xl2norm=sqrt(xl2norm*xl2norm+ul2*ul2);
            xnorm=sqrt(xl2norm*xl2norm+ul*ul+u*u);

            // Update w. Update x  except if it will become too big
            if(iter==1) {
                KERNEL_TYPE::copy(_wl,_wl2);        //wl2=wl;
                KERNEL_TYPE::copy(_v,_wl);
                KERNEL_TYPE::scale(sr1,_wl);        //wl=v*sr1;
                KERNEL_TYPE::copy(_v,_w);
                KERNEL_TYPE::scale(-cr1,_w);        //w=v*-cr1;
            } else if(iter==2) {
                KERNEL_TYPE::copy(_wl,_wl2);        //wl2=wl;
                KERNEL_TYPE::copy(_w,_wl);
                KERNEL_TYPE::scale(cr1,_wl);
                KERNEL_TYPE::addScaled(sr1,_v,_wl); //wl=w*cr1+v*sr1;// w=v_1
                KERNEL_TYPE::scale(sr1,_w);
                KERNEL_TYPE::addScaled(-cr1,_v,_w); //w=w*sr1+v*-cr1;
            } else {
                KERNEL_TYPE::copy(_wl,_wl2);        //wl2=wl;
                KERNEL_TYPE::copy(_w,_wl);          //wl=w;
                KERNEL_TYPE::copy(_wl2,_w);
                KERNEL_TYPE::scale(sr2,_w);
                KERNEL_TYPE::addScaled(-cr2,_v,_w); //w=wl2*sr2-v*cr2;
                KERNEL_TYPE::scale(cr2,_wl2);
                KERNEL_TYPE::addScaled(sr2,_v,_wl2);//wl2=wl2*cr2+v*sr2;
                KERNEL_TYPE::copy(_wl,_v);
                KERNEL_TYPE::scale(cr1,_v);
                KERNEL_TYPE::addScaled(sr1,_w,_v);  //v=wl*cr1+w*sr1;// temporarily using v for intermediate wl
                KERNEL_TYPE::scale(-cr1,_w);
                KERNEL_TYPE::addScaled(sr1,_wl,_w); //w=wl*sr1-w*cr1;
                KERNEL_TYPE::copy(_v,_wl);          //wl=v;
            }
            x1last=x[0];
            KERNEL_TYPE::addScaled(ul2,_wl2,_xl2);//xl2=xl2+wl2*ul2;// wl2=wl3;
            KERNEL_TYPE::copy(_xl2,x);
            KERNEL_TYPE::addScaled(ul,_wl,x);
            KERNEL_TYPE::addScaled(u,_w,x);//x=xl2+wl*ul+w*u;

            // Compute the next right plane rotation P{k-1,k+1}
            gamal_tmp=gamal;
            PMINRESSolver<T>::symGivens2(gamal,eplnn, cr2,sr2,gamal);

            // Estimate various norms
            rnorml=rnorm;   // ||r_{k-1}||
            Anorml=Anorm;
            Acondl=Acond;
            relresl=relres;
            L_ek_norml=L_ek_norm;

            if(flag != 9) {
                if(iter > 1) {
                    Tnorm=sqrt(Tnorm*Tnorm+beta*beta+alfa*alfa+betan*betan);   // Frobenius norm OF T_{k,k+1}
                    // used as estimate of AnormF.
                } else {
                    Tnorm=sqrt(alfa*alfa+betan*betan);
                }
                //AnormF=norm(Tnorm);
                L_el2_norm=sqrt(gamal2*gamal2+veplnl*veplnl+eta*eta);
                Anorm=std::max(Anorm,std::max(pnorm,gamal));
                rnorm=phi;                           // ||r_k||
                relres=rnorm/(Anorm*xnorm+beta1);    // ||r||/(||A||||x||+||b||)
                if(iter==1) {
                    gmin=gama;
                } else if(iter > 1) {
                    gmin=std::min(gmin,std::min(gamal,std::abs(gama)));
                    //gmin=min([gmin, gamal]);
                }
                Acond=Anorm/gmin;
            }

            rootl=sqrt(gbar*gbar+dltan*dltan);
            Arnorml=rnorml*rootl;                   // ||A r_{k-1} ||
            relAresl=rootl/Anorm;                   // ||Ar|| / (||A|| ||r||)
            //relAresl=Arnorml / Anorm;             // ||Ar|| / ||A||

            //---------------------------------------------------------------
            // See if any of the stopping criteria are satisfied.
            // In rare cases, flag is already -1 from above (Abar=const*I).
            //---------------------------------------------------------------
            //epsa=Anorm * eps;             % Used below to test betan and ||A r|| --- never used?
            epsx=Anorm*xnorm*eps;
            epsr=Anorm*xnorm*_toleranceFactor;
            // Test for singular Tk (hence singular A)
            // or x is already an LS solution (so again A must be singular).

            if((flag==0) || (flag==9)) {
                t1=1.0f+relres;                         // These tests work if rtol < eps
                t2=1.0f+relAresl;
                if(_toleranceFactor >=eps) {
                    if(t1 <= 1.0f)flag=3;                           // Accurate Ax=b solution
                    else if(t2 <= 1.0f)flag=4;                      // Accurate LS solution
                    else if(relres <= _toleranceFactor)flag=1;      // Good enough Ax=b solution
                    else if(relAresl <= _toleranceFactor)flag=2;    // Good enough LS solution
                    else if(epsx >= beta1)flag=5;                   // x is an eigenvector
                    else if(xnorm >= _maxxnorm)flag=6;              // xnorm exceeded its limit
                    else if(Acond >= _Acondlim)flag=7;              // Huge Acond (but maybe ok?)
                    else if(iter >= _maxIterations)flag=8;          // Too many itns
                } else {
                    if(relres <= _toleranceFactor)flag=1;           // Good enough Ax=b solution
                    else if(relAresl <= _toleranceFactor)flag=2;    // Good enough LS solution
                    else if(epsx >= beta1)flag=5;                   // x is an eigenvector
                    else if(xnorm >= _maxxnorm)flag=6;              // xnorm exceeded its limit
                    else if(Acond >= _Acondlim)flag=7;              // Huge Acond (but maybe ok?)
                    else if(iter >= _maxIterations)flag=8;          // Too many itns
                }
            }

            if(flag != 0) { // flag may have been changed above
                done=1;
                xnorm=sqrt(KERNEL_TYPE::dot(x,x));  // So far, xnorm is an estimate of sqrt(x'Mx)
                _fixedMatrix->multiply(x,_r1);
                KERNEL_TYPE::scale(-1.0f,_r1);
                KERNEL_TYPE::add(rhs,_r1,_r1);
                KERNEL_TYPE::addScaled(_shift,x,_r1);
                //r1 = b - minresxxxA( A, x , varargin ) + shift*x;

                relres=rnorm/(Anorm*xnorm+beta1);
                _fixedMatrix->multiply(_r1,_r2);
                KERNEL_TYPE::addScaled(-_shift,_r1,_r2);
                Arnorm=sqrt(KERNEL_TYPE::dot(_r2,_r2));
                //Arnorm=norm(minresxxxA( A, r1, varargin )-shift*r1);
                if(rnorm > numeric_limits<T>::min()) {
                    relAres=Arnorm/(Anorm*rnorm);
                }
                if((flag == 2) || (flag == 4) || (flag == 7)) {
                    // LS, or _maxxnorm exceeded, or Acondlim exceeded
                    lastiter=1;
                }
                if(lastiter || (flag == 6)) {
                    rnorm=rnorml;
                    relres=relresl;
                    if(lastiter) {
                        iter=iter-1;
                        pnorm=pnorml;
                    }
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
                         << std::setw(14) << relresl
                         << std::setw(14) << relAresl
                         << std::setw(14) << Anorml
                         << std::setw(14) << Acondl
                         << std::setw(14) << betal << std::endl;
            }
            _residualOut=relresl;
            if(_cb) {
                if((*_cb)(x,relresl,iter))
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
                     << std::setw(14) << relresl
                     << std::setw(14) << relAresl
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
protected:
    // internal structures
    PRECONDITIONER_TYPE _pre;
    Vec _xl2,_r1,_r2,_r3,_w,_wl,_wl2,_v; // temporary vectors for PCG
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
