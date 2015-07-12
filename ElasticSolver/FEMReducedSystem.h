#ifndef FEM_REDUCED_SYSTEM_H
#define FEM_REDUCED_SYSTEM_H

#include "FEMSystem.h"

PRJ_BEGIN

class ReducedMCalculator;
class ReducedKCalculator;
class FEMReducedSystem : public FEMSystem
{
public:
    enum KINETIC_MODEL_TYPE {
        //basis model reduction
        LMA,
        WARPPED_LMA,
        //simulate using RS coords
        POST_PROCESS_RS,
        COUPLED_RS,
        POST_PROCESS_REDUCED_RS,
        CUBATURE_COUPLED_REDUCED_RS,
        TAYLOR_COUPLED_REDUCED_RS,
        //verbose
        NR_REDUCED_SYSTEM_TYPE,
    };
    FEMReducedSystem(FEMBody& body);
    virtual boost::shared_ptr<FEMSystem> copy(FEMBody& other) const;
    virtual void setReducedModel(KINETIC_MODEL_TYPE type,sizeType nrCub,sizeType nrPose,scalar tolSC);
    //for integrator
    virtual void buildSystem(scalar MCoef,scalar KCoef,scalar CCoef,FEMSystemMatrix& LHS,FEMSystemMatrix& U,
                             const Vec& MRHSL,const Vec& CRHS,scalar FCoef,Vec& RHS,sizeType offVar) const;
    //IO
    virtual void getPos(Vec& P) const;
    virtual void getPosL(Vec& L) const;
    virtual Vec LtoF(const Vec& L) const;
    virtual void setPos(const Vec& P);
    virtual sizeType size() const;
    virtual sizeType sizeL() const;
    virtual void debugU();
    //(local) basis
    void writeBasisVTK(const std::string& path,scalar off,bool U=true);
    void buildU();
    void buildU(sizeType nr,bool overwrite=true);
    //user must make sure their current basis is LMA basis or the method fails
    virtual void onDirty();
protected:
    void removeRigidMode(SparseReducedBasis& basis,Vec& lambda) const;
    //reduced force calculation
    boost::shared_ptr<ReducedMCalculator> _MCalc;
    boost::shared_ptr<ReducedKCalculator> _KCalc;
};
class ReducedMCalculator
{
public:
    typedef Cold Vec;
    typedef vector<Eigen::Triplet<scalarD,sizeType> > TRIPS;
    ReducedMCalculator(FEMReducedSystem& system);
    virtual void buildBasis();
    virtual void precompute();
    virtual void init();
    void setPad(bool pad);
    //for integrator
    virtual void calculateM(const Vec& S,Matd& M,const Vec& MRHSL,Vec& RHS,sizeType r,sizeType c,FEMSystemMatrix* U);
    virtual Vec calculateFE(const Vec& S,Matd& K);
    //utility
    virtual void getPos(Vec& P) const;
    virtual void getPosL(Vec& L) const;
    virtual Vec LtoF(const Vec& L) const;
    virtual void setPos(const Vec& P);
    virtual sizeType size() const;
    virtual sizeType sizeL() const;
protected:
    FEMBody& _body;
    EnergyPool& _energys;
    //preprocess external energy terms
    Vec _UTf;
    Matd _UTHU;
};

PRJ_END

#endif