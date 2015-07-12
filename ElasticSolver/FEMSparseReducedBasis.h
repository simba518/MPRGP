#ifndef FEM_SPARSE_REDUCED_BASIS_H
#define FEM_SPARSE_REDUCED_BASIS_H

#include "FEMMesh.h"
#include <Eigen/Sparse>

PRJ_BEGIN

struct SparseReducedBasis : public Serializable {
    typedef Cold Vec;
    //constructor
    SparseReducedBasis();
    SparseReducedBasis(const SparseReducedBasis& other);
    SparseReducedBasis& operator=(const SparseReducedBasis& other);
    void resize(sizeType row,sizeType col);
    sizeType nrLMA(const Matd& B) const;
    //extended S
    boost::shared_ptr<Serializable> copy() const;
    void applyTrans(const Mat4& M);
    bool read(std::istream& is);
    bool write(std::ostream& os) const;
    //we have two basis U and B
    //basis U is for force-calculation
    //basis B is for coupledRS mass-calculation
    Matd _UTMU,_UTKU;
    Matd _U,_B;
    //high order polynomial coefficients
    Matd _STVKU,_RSP,_RSQ;
    Eigen::SparseMatrix<scalarD,0,sizeType> _M,_K;
};
struct Cubature : public Serializable {
    typedef Coli Veci;
    Cubature();
    boost::shared_ptr<Serializable> copy() const;
    bool read(std::istream& is);
    bool write(std::ostream& os) const;
    vector<FEMInterp> _tet;
    Cold _weight;
};

PRJ_END

#endif