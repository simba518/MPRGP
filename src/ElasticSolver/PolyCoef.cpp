#include "PolyCoef.h"
#include <boost/shared_ptr.hpp>

USE_PRJ_NAMESPACE

//polynomial coefficient extraction
PolyCoef::PolyCoef(Matd& coef,sizeType order):_coef(coef),_index(order) {}
void PolyCoef::findCoef(bool overwrite)
{
    //find index
    sizeType order=(sizeType)_index.size(),total=1;
    for(sizeType ord=0; ord<order; ord++) {
        _index[ord].clear();
        Coli tmp(ord+1);
        findIndex(_index[ord],tmp,0);
        total+=(sizeType)_index[ord].size();
    }
    //recompute
    if(overwrite || _coef.rows() != nrF() || _coef.cols() != total) {
        Matd RHS(total,nrF());
        Cold S(nrX()),f(nrF());
        TRIPS LHSTrips;
        LHSTrips.reserve(total*10);
        sizeType off=0;
        //add order 0
        S.setZero();
        addPatternLHS(LHSTrips,off,S);
        RHS.row(off++)=findF(S);
        //addPattern for order 1/2/3/4
        for(sizeType ord=0; ord<order; ord++) {
            sizeType nr=(sizeType)_index[ord].size();
            for(sizeType i=0; i<nr; i++,off++) {
                S.setZero();
                const Coli& index=_index[ord][i];
                for(sizeType j=0; j<index.size(); j++)
                    S[index[j]]+=1.0f;
                addPatternLHS(LHSTrips,off,S);
                RHS.row(off)=findF(S);
            }
        }
        ASSERT(off == total)
        //solve
        Eigen::SparseMatrix<scalarD,0,sizeType> LHS(total,total);
        LHS.setFromTriplets(LHSTrips.begin(),LHSTrips.end());
        Eigen::SparseLU<Eigen::SparseMatrix<scalarD,0,sizeType> > sol(LHS);
        ASSERT(sol.info() == Eigen::Success)
        _coef.resize(nrF(),total);
        for(sizeType i=0; i<nrF(); i++) {
            Cold tmp=sol.solve(RHS.col(i));
            _coef.row(i)=tmp;
        }
    }
}
Cold PolyCoef::findFK(const Cold& S,Matd& K) const
{
    sizeType order=(sizeType)_index.size(),total=0;
    //add order 0
    Cold F=_coef.col(total++);
    K.setZero(nrF(),nrX());
    //add order 1/2/3/4
    if(order >= 1)addOrder1(S,F,K,total);
    if(order >= 2)addOrder2(S,F,K,total);
    if(order >= 3)addOrder3(S,F,K,total);
    if(order >= 4)addOrder4(S,F,K,total);
    ASSERT_MSG(order <= 4,"We don't handle polynomial order > 4!")
    return F;
}
void PolyCoef::debugCoef() const
{
    for(sizeType i=0; i<10; i++) {
        Cold testS(nrX()),testS2;
        testS.setRandom();

        Matd KSub,KSub2;
        Cold fSub=findFK(testS,KSub),fSub2;
        INFOV("Func: %f %f",findF(testS).norm(),fSub.norm())
        ASSERT((findF(testS)-fSub).norm() < 1E-3f)

        for(sizeType j=0; j<nrX(); j++) {
#define DELTA 1E-7f
            testS2=testS+Cold::Unit(nrX(),j)*DELTA;
            fSub2=findFK(testS2,KSub2);
            INFOV("Jacobian: %f %f",((fSub2-fSub)/DELTA).norm(),KSub.col(j).norm());
            ASSERT(((fSub2-fSub)/DELTA-KSub.col(j)).norm() < 1E-3f)
#undef DELTA
        }
    }
}
void PolyCoef::findIndex(vector<Coli>& ind,Coli& curr,sizeType i) const
{
    if(i == curr.size()) {
        ind.push_back(curr);
        return;
    }
    sizeType from=(i==0)?0:curr[i-1];
    for(sizeType j=from; j<nrX(); j++) {
        curr[i]=j;
        findIndex(ind,curr,i+1);
    }
}
void PolyCoef::addOrder1(const Cold& S,Cold& F,Matd& K,sizeType& total) const
{
    sizeType nr=(sizeType)_index[0].size();
    Matd CBlk=_coef.block(0,total,nrF(),nr);
    Cold SBlk=S.block(0,0,nr,1);
    F+=CBlk*SBlk;
    K+=CBlk;
    total+=nr;
}
void PolyCoef::addOrder2(const Cold& S,Cold& F,Matd& K,sizeType& total) const
{
    Cold col;
    sizeType nr=(sizeType)_index[1].size();
    for(sizeType i=0,j=total; i<nr; i++,j++) {
        const Coli& index=_index[1][i];
        col=_coef.col(j);
        F+=col*(S[index[0]]*S[index[1]]);
        K.col(index[0])+=col*(S[index[1]]);
        K.col(index[1])+=col*(S[index[0]]);
    }
    total+=nr;
}
void PolyCoef::addOrder3(const Cold& S,Cold& F,Matd& K,sizeType& total) const
{
    Cold col;
    sizeType nr=(sizeType)_index[2].size();
    for(sizeType i=0,j=total; i<nr; i++,j++) {
        const Coli& index=_index[2][i];
        col=_coef.col(j);
        F+=col*(S[index[0]]*S[index[1]]*S[index[2]]);
        K.col(index[0])+=col*(S[index[1]]*S[index[2]]);
        K.col(index[1])+=col*(S[index[0]]*S[index[2]]);
        K.col(index[2])+=col*(S[index[0]]*S[index[1]]);
    }
    total+=nr;
}
void PolyCoef::addOrder4(const Cold& S,Cold& F,Matd& K,sizeType& total) const
{
    Cold col;
    sizeType nr=(sizeType)_index[3].size();
    for(sizeType i=0,j=total; i<nr; i++,j++) {
        const Coli& index=_index[3][i];
        col=_coef.col(j);
        F+=col*(S[index[0]]*S[index[1]]*S[index[2]]*S[index[3]]);
        K.col(index[0])+=col*(S[index[1]]*S[index[2]]*S[index[3]]);
        K.col(index[1])+=col*(S[index[0]]*S[index[2]]*S[index[3]]);
        K.col(index[2])+=col*(S[index[0]]*S[index[1]]*S[index[3]]);
        K.col(index[3])+=col*(S[index[0]]*S[index[1]]*S[index[2]]);
    }
    total+=nr;
}
void PolyCoef::addPatternLHS(TRIPS& trips,sizeType off,const Cold& S) const
{
    sizeType order=(sizeType)_index.size(),total=0;
    //add order 0
    trips.push_back(Eigen::Triplet<scalarD,sizeType>(off,total++,1.0f));
    //add order 1/2/3/4
    for(sizeType ord=0; ord<order; ord++) {
        sizeType nr=(sizeType)_index[ord].size();
        for(sizeType i=0; i<nr; i++,total++) {
            const Coli& index=_index[ord][i];
            scalarD val=1.0f;
            for(sizeType j=0; j<index.size(); j++)
                val*=S[index[j]];
            if(val>0.0f)
                trips.push_back(Eigen::Triplet<scalarD,sizeType>(off,total,val));
        }
    }
}

//precompute s
PolyCoefEfficient::PolyCoefEfficient(Matd& coef,sizeType order):_coef(coef),_order(order) {}
void PolyCoefEfficient::findCoef(bool overwrite)
{
    assemble();
    if(overwrite || _coef.rows() != nrF() || _coef.cols() != _off[4]) {
        Matd RHS(_off[4],nrF());
        Cold S(nrX()),f(nrF());
        TRIPS LHSTrips;
        LHSTrips.reserve(_off[4]*10);
        sizeType off=0;
        for(sizeType i=0; i<_off[4]; i++,off++) {
            S.setZero();
            for(sizeType t=0; t<_terms[i]._order; t++)
                S[_terms[i]._diffS[t]]+=1.0f;
            addPatternLHS(LHSTrips,off,S);
            RHS.row(off)=findF(S);
        }
        ASSERT(off==_off[4])
        //solve
        Eigen::SparseMatrix<scalarD,0,sizeType> LHS(_off[4],_off[4]);
        LHS.setFromTriplets(LHSTrips.begin(),LHSTrips.end());
        Eigen::SparseLU<Eigen::SparseMatrix<scalarD,0,sizeType> > sol(LHS);
        ASSERT(sol.info() == Eigen::Success)
        _coef.resize(nrF(),_off[4]);
        for(sizeType i=0; i<nrF(); i++) {
            Cold tmp=sol.solve(RHS.col(i));
            _coef.row(i)=tmp;
        }
    }
}
Cold PolyCoefEfficient::findFK(const Cold& S,Matd& K) const
{
    //Dynamic Programming
    sizeType nrTerm=(sizeType)_terms.size();
    vector<scalarD> val(nrTerm,1.0f);
    for(sizeType i=1; i<nrTerm; i++) {
        const Term& t=_terms[i];
        val[i]=(t._diffT[0]==-1)?S[t._diffS[0]]:val[t._diffT[0]]*S[t._diffS[0]];
    }
    //Assemble
    Cold ret=Cold::Zero(nrF());
    K.setZero(nrF(),nrX());
    Cold coeff;
    for(sizeType i=0; i<_off[4]; i++) {
        coeff=_coef.col(i);
        ret+=coeff*val[i];
        const Term& t=_terms[i];
        for(sizeType tt=0; tt<t._order; tt++)
            if(t._diffT[tt] == -1)
                K.col(t._diffS[tt])+=coeff;
            else K.col(t._diffS[tt])+=coeff*val[t._diffT[tt]];
    }
    return ret;
}
void PolyCoefEfficient::debugCoef() const
{
    for(sizeType i=0; i<10; i++) {
        Cold testS(nrX()),testS2;
        testS.setRandom();

        Matd KSub,KSub2;
        Cold fSub=findFK(testS,KSub),fSub2;
        INFOV("Func: %f %f",findF(testS).norm(),fSub.norm())
        ASSERT((findF(testS)-fSub).norm() < 1E-3f)

        for(sizeType j=0; j<nrX(); j++) {
#define DELTA 1E-7f
            testS2=testS+Cold::Unit(nrX(),j)*DELTA;
            fSub2=findFK(testS2,KSub2);
            INFOV("Jacobian: %f %f",((fSub2-fSub)/DELTA).norm(),KSub.col(j).norm());
            ASSERT(((fSub2-fSub)/DELTA-KSub.col(j)).norm() < 1E-3f)
#undef DELTA
        }
    }
}
void PolyCoefEfficient::assemble()
{
    //set term
    Term t;
    sizeType nr=nrX();
    sizeType nrTerm=0;
    _termMap.clear();
    _terms.clear();
    //zero order
    t._order=0;
    _terms.push_back(t);
    nrTerm++;
    //first order
    _off[0]=nrTerm;
    if(_order >= 1)
        for(sizeType i=0; i<nr; i++,nrTerm++) {
            _termMap[getKey(Vec4i(i,-1,-1,-1))]=nrTerm;
            t._order=1;
            t._diffS[0]=i;
            t._diffT[0]=-1;
            _terms.push_back(t);
        }
    //second order
    _off[1]=nrTerm;
    if(_order >= 2)
        for(sizeType i=0; i<nr; i++)
            for(sizeType j=i; j<nr; j++,nrTerm++) {
                _termMap[getKey(Vec4i(i,j,-1,-1))]=nrTerm;
                t._order=2;
                t._diffS[0]=j;
                t._diffT[0]=_termMap.find(getKey(Vec4i(i,-1,-1,-1)))->second;
                t._diffS[1]=i;
                t._diffT[1]=_termMap.find(getKey(Vec4i(j,-1,-1,-1)))->second;
                _terms.push_back(t);
            }
    //third order
    _off[2]=nrTerm;
    if(_order >= 3)
        for(sizeType i=0; i<nr; i++)
            for(sizeType j=i; j<nr; j++)
                for(sizeType k=j; k<nr; k++,nrTerm++) {
                    _termMap[getKey(Vec4i(i,j,k,-1))]=nrTerm;
                    t._order=3;
                    t._diffS[0]=k;
                    t._diffT[0]=_termMap.find(getKey(Vec4i(i,j,-1,-1)))->second;
                    t._diffS[1]=j;
                    t._diffT[1]=_termMap.find(getKey(Vec4i(i,k,-1,-1)))->second;
                    t._diffS[2]=i;
                    t._diffT[2]=_termMap.find(getKey(Vec4i(j,k,-1,-1)))->second;
                    _terms.push_back(t);
                }
    //forth order
    _off[3]=nrTerm;
    if(_order >= 4)
        for(sizeType i=0; i<nr; i++)
            for(sizeType j=i; j<nr; j++)
                for(sizeType k=j; k<nr; k++)
                    for(sizeType l=k; l<nr; l++,nrTerm++) {
                        _termMap[getKey(Vec4i(i,j,k,l))]=nrTerm;
                        t._order=4;
                        t._diffS[0]=l;
                        t._diffT[0]=_termMap.find(getKey(Vec4i(i,j,k,-1)))->second;
                        t._diffS[1]=k;
                        t._diffT[1]=_termMap.find(getKey(Vec4i(i,j,l,-1)))->second;
                        t._diffS[2]=j;
                        t._diffT[2]=_termMap.find(getKey(Vec4i(i,k,l,-1)))->second;
                        t._diffS[3]=i;
                        t._diffT[3]=_termMap.find(getKey(Vec4i(j,k,l,-1)))->second;
                        _terms.push_back(t);
                    }
    _off[4]=nrTerm;
}
sizeType PolyCoefEfficient::getKey(const Vec4i& index) const
{
    unsigned short val[4];
    for(char i=0; i<4; i++)
        val[i]=((unsigned short)index[i]==-1) ? 65535 : (unsigned short)index[i];
    return *((sizeType*)val);
}
Vec4i PolyCoefEfficient::setKey(sizeType index) const
{
    Vec4i ret;
    unsigned short *val=(unsigned short*)&index;
    for(char i=0; i<4; i++)
        ret[i]=val[i] == 65535 ? -1 : val[i];
    return ret;
}
void PolyCoefEfficient::addPatternLHS(TRIPS& trips,sizeType off,const Cold& S) const
{
    for(sizeType i=0; i<_off[4]; i++) {
        const Term& t=_terms[i];
        scalarD val=1.0f;
        for(sizeType tt=0; tt<t._order && val>0.0f; tt++)
            val*=S[t._diffS[tt]];
        if(val>0.0f)
            trips.push_back(Eigen::Triplet<scalarD,sizeType>(off,i,val));
    }
}

//debugger
struct DebugPoly {
public:
    DebugPoly(sizeType nr) {
        _H1.setRandom(nr,nr);
        _H2.setRandom(nr,nr);
        _A1.setRandom(nr);
        _A2.setRandom(nr);
        _b=rand()/(scalarD)RAND_MAX;
    }
    virtual scalarD findF(const Cold& S) const {
        return (_A1.dot(S)+S.dot(_H1*S)*0.5f)*(_A2.dot(S)+S.dot(_H2*S)*0.5f)+_b;
    }
    Matd _H1,_H2;
    Cold _A1,_A2;
    scalarD _b;
};
template <typename PARENT,int NP>
struct DebugNPoly : public PARENT {
    DebugNPoly(Matd& coef,sizeType nr):PARENT(coef,4) {
        _nPoly.resize(NP);
        for(sizeType i=0; i<NP; i++)
            _nPoly[i].reset(new DebugPoly(nr));
    }
    virtual Cold findF(const Cold& S) const {
        Cold ret(nrF());
        for(sizeType i=0; i<nrF(); i++)
            ret[i]=_nPoly[i]->findF(S);
        return ret;
    }
    virtual sizeType nrF() const {
        return NP;
    }
    virtual sizeType nrX() const {
        return _nPoly[0]->_A1.size();
    }
    vector<boost::shared_ptr<DebugPoly> > _nPoly;
};
void PolyCoefEfficient::debug()
{
    Matd coef;
    DebugNPoly<PolyCoefEfficient,100> dbg(coef,10);
    dbg.findCoef(false);
    dbg.debugCoef();
}
void PolyCoef::debug()
{
    Matd coef;
    DebugNPoly<PolyCoef,100> dbg(coef,10);
    dbg.findCoef(false);
    dbg.debugCoef();
}