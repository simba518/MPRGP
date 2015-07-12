#ifndef FEM_ENERGY_H
#define FEM_ENERGY_H

#include "MathBasic.h"
#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>
#include <boost/property_tree/ptree.hpp>
#include <Eigen/Sparse>
#include <typeinfo>

PRJ_BEGIN

struct FEMVertex;
struct FEMCell;
struct FEMInterp;
struct FEMBody;
class FEMMesh;

struct BlockSparseMatrix;
struct BasisExtractor {
    virtual Matd operator()(sizeType i) const=0;
    virtual Vec3d operator[](sizeType i) const=0;
};
struct MatrixBasisExtractor : public BasisExtractor {
    MatrixBasisExtractor(const Matd& U,const Cold& S);
    virtual Matd operator()(sizeType i) const;
    virtual Vec3d operator[](sizeType i) const;
    const Matd& _U;
    Cold _S;
};
class FEMEnergy
{
public:
    typedef Cold Vec;
    typedef vector<Eigen::Triplet<scalarD,sizeType> > TRIPS;
    virtual scalarD eval(Vec* f,TRIPS* H,scalar CF,scalar CH) const=0;
    virtual scalarD eval(const BasisExtractor& U,Vec* f,Matd* H,scalar CF,scalar CH) const;
    virtual void flagVertex(vector<bool>& flag) const;
    virtual boost::shared_ptr<FEMEnergy> copy(const FEMBody& other) const;
    void debugEnergy(const FEMBody& body);
protected:
    static void evalHelper(const BasisExtractor& U,Vec* f,Matd* H,Eigen::Matrix<scalarD,12,1>* fB,Eigen::Matrix<scalarD,12,12>* HB,const FEMCell& cell);
    static void evalHelper(Vec* f,TRIPS* H,Eigen::Matrix<scalarD,12,1>* fB,Eigen::Matrix<scalarD,12,12>* HB,const FEMCell& cell);
};
class FixedForceEnergy : public FEMEnergy
{
public:
    FixedForceEnergy(boost::shared_ptr<FEMVertex> v,Vec3 a);
    virtual scalarD eval(Vec* f,TRIPS* H,scalar CF,scalar CH) const;
    virtual scalarD eval(const BasisExtractor& U,Vec* f,Matd* H,scalar CF,scalar CH) const;
    virtual void flagVertex(vector<bool>& flag) const;
    virtual boost::shared_ptr<FEMEnergy> copy(const FEMBody& other) const;
protected:
    boost::shared_ptr<FEMVertex> _v;
    Vec3 _f;
};
class FixVertexEnergy : public FEMEnergy
{
public:
    FixVertexEnergy(sizeType dim,boost::shared_ptr<FEMVertex> v,scalar K,const Vec3& p0);
    virtual scalarD eval(Vec* f,TRIPS* H,scalar CF,scalar CH) const;
    virtual scalarD eval(const BasisExtractor& U,Vec* f,Matd* H,scalar CF,scalar CH) const;
    virtual void flagVertex(vector<bool>& flag) const;
    virtual boost::shared_ptr<FEMEnergy> copy(const FEMBody& other) const;
    Vec3 getPos() const;
    void setPos(const Vec3& pos);
    boost::shared_ptr<FEMVertex> getVert() const;
protected:
    //data
    sizeType _dim;
    Vec3 _pos;
    scalar _K;
    //geometry
    boost::shared_ptr<FEMVertex> _v;
};
class FixPointEnergy : public FEMEnergy
{
public:
    FixPointEnergy(sizeType dim,boost::shared_ptr<FEMCell> c,const FEMInterp& I,scalar K,const Vec3& p0);
    virtual scalarD eval(Vec* f,TRIPS* H,scalar CF,scalar CH) const;
    virtual scalarD eval(const BasisExtractor& U,Vec* f,Matd* H,scalar CF,scalar CH) const;
    virtual void flagVertex(vector<bool>& flag) const;
    virtual boost::shared_ptr<FEMEnergy> copy(const FEMBody& other) const;
    Vec3 getPos() const;
    void setPos(const Vec3& pos);
    const FEMInterp& getI() const;
    boost::shared_ptr<FEMCell> getCell() const;
protected:
    //data
    sizeType _dim;
    Vec3 _pos;
    scalar _K;
    Eigen::Matrix<scalarD,12,12> _HB;
    //geometry
    boost::shared_ptr<FEMInterp> _I;
    boost::shared_ptr<FEMCell> _c;
};
class MaterialEnergy : public FEMEnergy
{
public:
    enum TYPE {
        STVK,
        LINEAR,
        COROTATIONAL,
        COROTATIONAL_EXACT,
        NONHK,
        FUNG,
    };
    MaterialEnergy();
    MaterialEnergy(boost::shared_ptr<FEMCell> c,const boost::property_tree::ptree& tree);
    virtual scalarD eval(Vec* f,TRIPS* H,scalar CF,scalar CH) const;
    virtual scalarD eval(const BasisExtractor& U,Vec* f,Matd* H,scalar CF,scalar CH) const;
    virtual scalarD eval(Eigen::Matrix<scalarD,12,1>& fB,Eigen::Matrix<scalarD,12,12>* HB) const;
    virtual void flagVertex(vector<bool>& flag) const;
    boost::shared_ptr<FEMEnergy> copy(const FEMBody& other) const;
    void setInvertible(bool invertible);
    TYPE getType() const;
    void setType(TYPE type);
    scalar getVol() const;
    //for isotropic-model
    scalar getMu() const;
    void setMu(scalar mu);
    scalar getLambda() const;
    void setLambda(scalar lambda);
    //for Fung's soft-tissue model
    scalar getMuFung() const;
    void setMuFung(scalar muFung);
    scalar getLambdaFung() const;
    void setLambdaFung(scalar lambdaFung);
    scalar getCFung() const;
    void setCFung(scalar cFung);
    boost::shared_ptr<FEMCell> getCell() const;
protected:
    void setParam(const boost::property_tree::ptree& tree);
    static scalarD evalF(const Mat3& f,scalarD lambda,scalarD mu,scalarD V,scalarD dimMask,
                         Eigen::Matrix<scalarD,9,1>& grad,Eigen::Matrix<scalarD,9,9>* hess);
    static scalarD evalFLinear(const Mat3& f,scalarD lambda,scalarD mu,scalarD V,scalarD dimMask,
                               Eigen::Matrix<scalarD,9,1>& grad,Eigen::Matrix<scalarD,9,9>* hess);
    static scalarD evalFCorot(const Mat3& f,const Mat3& fR,
                              scalarD lambda,scalarD mu,scalarD V,scalarD dimMask,
                              Eigen::Matrix<scalarD,9,1>& grad,Eigen::Matrix<scalarD,9,9>* hess);
    static scalarD evalFCorotE(const Mat3& f,const Mat3& fR,const Eigen::Matrix<scalarD,9,9>& DRDF,
                               scalarD lambda,scalarD mu,scalarD V,scalarD dimMask,
                               Eigen::Matrix<scalarD,9,1>& grad,Eigen::Matrix<scalarD,9,9>* hess);
    static scalarD evalFCorotE(const Mat3& f,const Mat3& fR,scalarD lambda,scalarD mu,scalarD V,scalarD dimMask,
                               Eigen::Matrix<scalarD,9,1>& grad,Eigen::Matrix<scalarD,9,18>* hess);
    static scalarD evalFNonHK(const Mat3& f,scalarD lambda,scalarD mu,scalarD V,scalarD dimMask,scalarD minVol,
                              Eigen::Matrix<scalarD,9,1>& grad,Eigen::Matrix<scalarD,9,9>* hess);
    static scalarD evalFFung(const Mat3& f,scalarD lambda,scalarD mu,scalarD lambda2,scalarD mu2,scalarD c,scalarD V,scalarD dimMask,
                             Eigen::Matrix<scalarD,9,1>& grad,Eigen::Matrix<scalarD,9,9>* hess);
    boost::shared_ptr<FEMCell> _c;
    scalar _mu,_lambda,_dimMask;
    scalar _muFung,_lambdaFung,_cFung;
    Eigen::Matrix<scalarD,9,12> _GComp;
    TYPE _type;
    bool _invertible;
};
class EnergyPool
{
public:
    typedef Cold Vec;
    typedef vector<Eigen::Triplet<scalarD,sizeType> > TRIPS;
    typedef vector<boost::shared_ptr<FEMEnergy> > ELIST;
    struct structE {
        structE():_lock(0) {}
        sizeType _lock;
        ELIST _elist;
    };
    typedef boost::unordered_map<std::string,structE> EMAP;
    typedef boost::unordered_map<std::string,sizeType> SMAP;
    typedef vector<boost::shared_ptr<FEMEnergy> >::const_iterator ConstIterator;
    struct TypeChooser {
        virtual bool operator()(const type_info& t) const {
            return true;
        }
    };
    struct NonMatTypeChooser : public TypeChooser {
        virtual bool operator()(const type_info& t) const {
            return t != typeid(MaterialEnergy);
        }
    };
    struct MatTypeChooser : public TypeChooser {
        virtual bool operator()(const type_info& t) const {
            return t == typeid(MaterialEnergy);
        }
    };
    struct FixPointEnergyChooser : public TypeChooser {
        virtual bool operator()(const type_info& t) const {
            return t == typeid(FixPointEnergy) || t == typeid(FixVertexEnergy);
        }
    };
public:
    EnergyPool(const FEMBody& body);
    ConstIterator begin(const type_info& t) const;
    ConstIterator end(const type_info& t) const;
    void addEnergy(boost::shared_ptr<FEMEnergy> e);
    void delEnergy(boost::shared_ptr<FEMEnergy> e);
    void delEnergy(const std::string& name,sizeType id);
    void setMaterialType(MaterialEnergy::TYPE type);
    MaterialEnergy::TYPE getMaterialType() const;
    void lock(bool setLock,const TypeChooser& c=TypeChooser());
    void lockNonMaterial(bool setLock) {
        lock(setLock,NonMatTypeChooser());
    }
    void clear();
    void clear(const std::string& name);
    void copy(const EnergyPool& other);
    scalarD energyEval(Vec* f,TRIPS* H,scalar CF,scalar CH,const TypeChooser& c=TypeChooser(),bool lock=true) const;
    scalarD energyEval(const BasisExtractor& U,Vec* f,Matd* H,scalar CF,scalar CH,const TypeChooser& c=TypeChooser(),bool lock=true) const;
    scalarD nonMaterialEnergyEval(Vec* f,TRIPS* H,scalar CF,scalar CH,bool lock=true) const {
        return energyEval(f,H,CF,CH,NonMatTypeChooser(),lock);
    }
    scalarD nonMaterialEnergyEval(const BasisExtractor& U,Vec* f,Matd* H,scalar CF,scalar CH,bool lock=true) const {
        return energyEval(U,f,H,CF,CH,NonMatTypeChooser(),lock);
    }
    scalarD materialEnergyEval(Vec* f,TRIPS* H,scalar CF,scalar CH,bool lock=true) const {
        return energyEval(f,H,CF,CH,MatTypeChooser(),lock);
    }
    scalarD materialEnergyEval(const BasisExtractor& U,Vec* f,Matd* H,scalar CF,scalar CH,bool lock=true) const {
        return energyEval(U,f,H,CF,CH,MatTypeChooser(),lock);
    }
    void flagVertex(vector<bool>& flag,const TypeChooser& c=TypeChooser()) const;
    void debugEnergy();
private:
    const FEMBody& _body;
    EMAP _energys;
};

PRJ_END

#endif
