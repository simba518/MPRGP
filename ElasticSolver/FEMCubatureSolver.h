#ifndef FEM_CUBATURE_SOLVER_H
#define FEM_CUBATURE_SOLVER_H

#include "FEMMesh.h"
#include "FEMSparseReducedBasis.h"
#include "solvers/MatVec.h"

PRJ_BEGIN

//FEM cubature solver solves the following problem:
//for an FEMMesh
//for a set N poses {s_i} in reduced basis U
//for an nonlinear function F(s,x)
//(where x is the rest pose tetrahedron vertex position)
//you need to find:
//|\Sigma_c w_c F(s_i,x_c) - \Sigma_t F(s_i,x_t)|_2^2 + \phi(w_c)
//(where \phi is a sparse regularizer)
//we use the GreedySparseCodingSolver to find a set of initial tets
//we run an EM algorithm of dictionary learning to jointly optimize {w_c,x_c}
class FEMCubatureProb
{
public:
    typedef Cold Vec;
    virtual ~FEMCubatureProb() {}
    virtual Matd operator()(const sizeType cid,const Matd& S) const=0;
    virtual void debugCallback(const Vec& b,const Matd& S) const {}
    virtual Vec2i nr() const=0;
    //Utility function: user can transform the metric of M according to W in fullspace defined by B
    //to do this, we find a B1 such that span(B1)=span(B) and B1^T*W*B1=I
    //and we solve argmin_s|B1*s-B*M*s0|, this amount to set: s=(B1^T*B1)^{-1}B1^T*B*M*s0
    //in the end, we set: M=(B1^T*B1)^{-1}B1^T*B*M
    static void transformMetric(const Eigen::SparseMatrix<scalarD,0,sizeType>& W,const Matd& B,Matd& M);
};
class FEMCubatureSolver
{
public:
    typedef Cold Vec;
    FEMCubatureSolver(FEMBody& body,FEMCubatureProb& obj);
    virtual ~FEMCubatureSolver() {}
    void generate(const Vec& weight,bool debugWrite=true);
    void solve(boost::shared_ptr<Cubature>& ptrCub,const std::string& extension);
    void solveSparseCodingHier(vector<std::pair<Coli,Vec> >& ret);
    void solveSparseCoding(Matd& A,vector<std::pair<Coli,Vec> >& ret,sizeType nr,vector<sizeType>* keyMap) const;
    const Matd& getS() const;
protected:
    bool read(std::istream& is);
    void scaleCol(std::iostream& os) const;
    void groupCol(Matd& AGroup,std::istream& is,const vector<sizeType>& keyMap,sizeType nrG) const;
    void selectCol(Matd& AGroup,std::istream& is,const vector<sizeType>* keyMap) const;
    void debugOutOfCore();
    //data
    FEMBody& _body;
    FEMCubatureProb& _obj;
    //problem
    Vec2i _szA;
    Vec _b,_w,_coef;
    Matd _S;
};
struct CFEMBody : public FEMBody {
    friend class FEMCubatureSolver;
    CFEMBody();
    CFEMBody(const FEMBody& body,const Cubature& ce);
    virtual FEMBody& operator=(const FEMBody& other);
    boost::shared_ptr<Serializable> copy() const;
    Eigen::SparseMatrix<scalarD,0,sizeType> _cInterp;
};

PRJ_END

#endif