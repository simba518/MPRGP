#ifndef FEM_SYSTEM_H
#define FEM_SYSTEM_H

#include "FEMEnergy.h"
#include "FEMSparseReducedBasis.h"
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

PRJ_BEGIN

struct ReducedBasis;
struct FEMBody;
class FEMMesh;
class FEMGeom;
class FEMFrameData
{
public:
    typedef Cold Vec;
    bool read(std::istream& is);
    bool write(std::ostream& os) const;
    Vec _X,_VL,_AL;
};
class FEMSystemMatrix
{
public:
    typedef Cold Vec;
    typedef vector<Eigen::Triplet<scalarD,sizeType> > TRIPS;
    FEMSystemMatrix();
    bool isDense() const;
    Vec solve(const Vec& RHS) const;
    void getDense(Matd& dense) const;
  
    int rows()const{
	  return _size[0];
    }
    int cols()const{
	  return _size[1];
    }
    const boost::shared_ptr<const TRIPS> getSparse()const{
	  return _sparse;
	}
    
    void getSparse(Eigen::SparseMatrix<scalarD,0,sizeType>& sparse) const;
    void clear();
    void reset(sizeType r,sizeType c,bool dense);
    void addBlock(sizeType r,sizeType c,const Matd& block);
    void addBlock(sizeType r,sizeType c,const TRIPS& trips,scalarD coef);
    void addI(sizeType r,sizeType c,const Vec& diag,scalarD coef);
    void addDense(sizeType r,sizeType c,const Matd& other);
    void addSparse(sizeType r,sizeType c,const Eigen::SparseMatrix<scalarD,0,sizeType>& other);
protected:
    //data
    boost::shared_ptr<TRIPS> _sparse;
    boost::shared_ptr<Matd> _dense;
    Vec2i _size;
};
class FEMSystem
{
public:
    typedef Cold Vec;
    typedef vector<Eigen::Triplet<scalarD,sizeType> > TRIPS;
    //constructor
    FEMSystem(FEMBody& body);
    virtual ~FEMSystem() {}
    FEMBody& body();
    virtual void setDamping(scalar alpha,scalar beta);
    virtual void resetParam(scalar alpha,scalar beta,scalar rho);
    virtual boost::shared_ptr<FEMSystem> copy(FEMBody& other) const;
    virtual void applyTrans(const Mat4& R);
    virtual void beforeCollision();
    virtual void onCollision(sizeType id);
    virtual void onCollisionSelf(sizeType id);
    virtual void afterCollision();
    //for integrator
    virtual void getSystemEnergy(scalarD& EK,scalarD& EP) const;
    virtual void buildSystem(scalar MCoef,scalar KCoef,scalar CCoef,FEMSystemMatrix& LHS,FEMSystemMatrix& U,
                             const Vec& MRHSL,const Vec& CRHS,scalar FCoef,Vec& RHS,sizeType offVar) const;
    virtual void buildU(Eigen::SparseMatrix<scalarD,0,sizeType>& U);
    virtual sizeType size() const;
    virtual sizeType sizeL() const;
    bool selfColl() const;
    virtual void debugU();
    //IO
    virtual void getFrame(FEMFrameData& data);
    virtual void setFrame(const FEMFrameData& data);
    virtual void getPos(Vec& P) const;
    virtual void getPosL(Vec& L) const;
    virtual void getVelL(Vec& VL) const;
    virtual void getAccelL(Vec& AL) const;
    virtual Vec LtoF(const Vec& F) const;
    virtual void setPos(const Vec& P);
    virtual void setVelL(const Vec& VL);
    virtual void setAccelL(const Vec& AL);
    //energy utility
    void addEnergyMaterial(scalar young,scalar poisson,MaterialEnergy::TYPE type,bool invertible=true);
    void setMaterialType(MaterialEnergy::TYPE type);
    void addEnergyMass();
    void addEnergyMass(const Vec3& g,const BBox<scalar>* bbm=NULL);
    void addEnergy(boost::shared_ptr<FEMEnergy> e);
    void delEnergy(boost::shared_ptr<FEMEnergy> e);
    void addForcePoint(sizeType id,const Vec3& force);
    void delForcePoint(sizeType id);
    EnergyPool& energys();
    void clearEnergy();
    void clearEnergy(const std::string& name);
    void readEnergy(const std::string& path);
    void readEnergy(const std::string& path,MaterialEnergy::TYPE type,bool invertible=true,bool constraints=true);
    void writeEnergy(const std::string& pathIn,const std::string& pathOut) const;
    //constraint utility
    void setConstraints(const boost::unordered_map<sizeType,Vec3>& cons);
    void addConstraintPoint(FEMInterp p,scalar K);
    void addConstraintPoint(sizeType id, const Vec3& pos);
    void addConstraintPoint(sizeType id);
    void delConstraintPoint(sizeType id);
    boost::unordered_map<sizeType,Vec3>& constraints();
    void clearConstraint();
    virtual void onDirty();
    virtual bool removeConstraint(TRIPS* LHS,Vec* RHS,boost::unordered_set<sizeType>* fixSetRet,sizeType offsetLHS) const;
protected:
    void applyConstraint();
    void clearVel(sizeType id);
    void clearAccel(sizeType id);
    //data
    FEMBody& _body;
    Vec _rho,_M;
    Vec _VL,_AL;
    boost::unordered_map<sizeType,Vec3> _constraints;
    EnergyPool _energys;
    bool _dirty,_selfColl;
};
class FEMSolver
{
public:
    typedef Cold Vec;
    typedef vector<Eigen::Triplet<scalarD,sizeType> > TRIPS;
    //constructor
    FEMSolver(sizeType dim,sizeType cOption=1);
    template <typename T>T& getSystem(sizeType b);
    void resetParam(scalar beta,scalar beta2,scalar gamma,scalar eps=1E-4f,sizeType maxIter=10);
    void resetImplicitEuler(scalar eps=1E-4f,sizeType maxIter=10) {
        resetParam(1.0f,0.5f,1.0f,eps,maxIter);
    }
    void resetCentralDifference(scalar eps=1E-4f,sizeType maxIter=10) {
        resetParam(0.5f,0.5f,0.5f,eps,maxIter);
    }
    void setCollK(scalar collK);
    FEMMesh& getMesh();
    const FEMMesh& getMesh() const;
    void setSelfColl(bool selfColl);
    //solver utility
    void advance(scalar dt);
    void resetBody(sizeType bid);
    void getSystemEnergy(scalarD& EK,scalarD& EP) const;
    //IO
    void writeFrame(std::ostream& os);
    void readFrame(std::istream& is);
    //energy utility
    void addConstraintPoint(FEMGeom& geom);
    void addConstraintPoint(const BBox<scalar>& bb);
    void addEnergy(boost::shared_ptr<FixPointEnergy> e);
    void delEnergy(boost::shared_ptr<FixPointEnergy> e);
    //geometry
    boost::shared_ptr<FEMGeom> _geom;
    boost::property_tree::ptree _tree;
protected:
    void assembleOffVar();
    sizeType nrVar() const;
    sizeType nrVarL() const;
    sizeType nrVarF() const;
    Vec LtoF(const Vec& L) const;
    boost::shared_ptr<FEMMesh> _mesh;
    vector<sizeType> _offVar,_offVarL,_offVarF;
    FEMSystemMatrix _LHS,_U;
};

PRJ_END

#include "FEMMesh.h"
PRJ_BEGIN
template <typename T>T& FEMSolver::getSystem(sizeType b)
{
    return boost::dynamic_pointer_cast<T>(_mesh->getB(b)._system).operator*();
}
PRJ_END

#endif
