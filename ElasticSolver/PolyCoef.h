#ifndef POLY_COEF_H
#define POLY_COEF_H

#include "MathBasic.h"
#include <Eigen/Sparse>
#include <boost/unordered_map.hpp>

PRJ_BEGIN

//polynomial coefficient extraction
struct PolyCoef {
    typedef vector<Eigen::Triplet<scalarD,sizeType> > TRIPS;
    PolyCoef(Matd& coef,sizeType order);
    void findCoef(bool overwrite);
    Cold findFK(const Cold& S,Matd& K) const;
    void debugCoef() const;
    static void debug();
    virtual Cold findF(const Cold& S) const=0;
    virtual sizeType nrF() const=0;
    virtual sizeType nrX() const=0;
protected:
    void findIndex(vector<Coli>& ind,Coli& curr,sizeType i) const;
    void addOrder1(const Cold& S,Cold& F,Matd& K,sizeType& total) const;
    void addOrder2(const Cold& S,Cold& F,Matd& K,sizeType& total) const;
    void addOrder3(const Cold& S,Cold& F,Matd& K,sizeType& total) const;
    void addOrder4(const Cold& S,Cold& F,Matd& K,sizeType& total) const;
    void addPatternLHS(TRIPS& trips,sizeType off,const Cold& S) const;
    //data
    vector<vector<Coli> > _index;
    Matd& _coef;
};
//more efficient version
struct PolyCoefEfficient {
    typedef vector<Eigen::Triplet<scalarD,sizeType> > TRIPS;
    struct Term {
        sizeType _order;
        sizeType _diffS[4],_diffT[4];
    };
    PolyCoefEfficient(Matd& coef,sizeType order);
    void findCoef(bool overwrite);
    Cold findFK(const Cold& S,Matd& K) const;
    void debugCoef() const;
    static void debug();
    virtual Cold findF(const Cold& S) const=0;
    virtual sizeType nrF() const=0;
    virtual sizeType nrX() const=0;
protected:
    void assemble();
    void addPatternLHS(TRIPS& trips,sizeType off,const Cold& S) const;
    sizeType getKey(const Vec4i& index) const;
    Vec4i setKey(sizeType index) const;
    boost::unordered_map<sizeType,sizeType> _termMap;
    vector<Term> _terms;
    sizeType _off[5],_order;
    Matd& _coef;
};

PRJ_END

#endif