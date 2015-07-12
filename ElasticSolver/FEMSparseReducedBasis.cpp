#include "FEMSparseReducedBasis.h"
#include "solvers/MatVec.h"

USE_PRJ_NAMESPACE

SparseReducedBasis::SparseReducedBasis():Serializable(5) {}
SparseReducedBasis::SparseReducedBasis(const SparseReducedBasis& other):Serializable(5)
{
    operator=(other);
}
SparseReducedBasis& SparseReducedBasis::operator=(const SparseReducedBasis& other)
{
    //kinetic-potential
    _UTMU=other._UTMU;
    _UTKU=other._UTKU;
    _U=other._U;
    _B=other._B;
    //tensor
    _STVKU=other._STVKU;
    _RSP=other._RSP;
    _RSQ=other._RSQ;
    return *this;
}
void SparseReducedBasis::resize(sizeType row,sizeType col)
{
    //kinetic-potential
    _UTMU.resize(col,col);
    _UTMU.setZero();
    _UTKU.resize(col,col);
    _UTKU.setZero();
    _U.setZero(row,col);
    _B.setZero(row,col);
    //matrix
    _M.resize(row,row);
    _K.resize(row,row);
}
sizeType SparseReducedBasis::nrLMA(const Matd& B) const
{
    Matd UTKU=B.transpose()*(_K*B).eval();
    sizeType nrLMA=0;
    scalarD thres=1E-4f;
    for(sizeType i=0; i<UTKU.cols(); i++,nrLMA++) {
        Vec col=UTKU.block(0,i,i+1,1);
        scalarD diag=col[i];
        col[i]=0.0f;
        if(col.norm()/diag > thres)
            break;
    }
    return nrLMA;
}
boost::shared_ptr<Serializable> SparseReducedBasis::copy() const
{
    return boost::shared_ptr<Serializable>(new SparseReducedBasis);
}
void SparseReducedBasis::applyTrans(const Mat4& M)
{
    for(sizeType i=0; i<_U.rows(); i+=3)
        _U.block(i,0,3,_U.cols())=
            M.block<3,3>(0,0).cast<scalarD>()*_U.block(i,0,3,_U.cols());
    for(sizeType i=0; i<_B.rows(); i+=3)
        _B.block(i,0,3,_B.cols())=
            M.block<3,3>(0,0).cast<scalarD>()*_B.block(i,0,3,_B.cols());
}
bool SparseReducedBasis::read(std::istream& is)
{
    //basis
    //kinetic-potential
    readBinaryData(_UTMU,is);
    readBinaryData(_UTKU,is);
    readBinaryData(_U,is);
    readBinaryData(_B,is);
    //tensor
    readBinaryData(_STVKU,is);
    readBinaryData(_RSP,is);
    readBinaryData(_RSQ,is);
    FixedSparseMatrix<scalarD,Kernel<scalarD> > K;
    K.read(is);
    K.toEigen(_K);
    FixedSparseMatrix<scalarD,Kernel<scalarD> > M;
    M.read(is);
    M.toEigen(_M);
    return is.good();
}
bool SparseReducedBasis::write(std::ostream& os) const
{
    //kinetic-potential
    writeBinaryData(_UTMU,os);
    writeBinaryData(_UTKU,os);
    writeBinaryData(_U,os);
    writeBinaryData(_B,os);
    //tensor
    writeBinaryData(_STVKU,os);
    writeBinaryData(_RSP,os);
    writeBinaryData(_RSQ,os);
    FixedSparseMatrix<scalarD,Kernel<scalarD> > K;
    K.fromEigen(_K);
    K.write(os);
    FixedSparseMatrix<scalarD,Kernel<scalarD> > M;
    M.fromEigen(_M);
    M.write(os);
    return os.good();
}

Cubature::Cubature():Serializable(6) {}
boost::shared_ptr<Serializable> Cubature::copy() const
{
    return boost::shared_ptr<Serializable>(new Cubature);
}
bool Cubature::read(std::istream& is)
{
    IOData dat;
    dat.registerType<FEMInterp>();
    readVector(_tet,is,&dat);
    readBinaryData(_weight,is);
    return is.good();
}
bool Cubature::write(std::ostream& os) const
{
    IOData dat;
    dat.registerType<FEMInterp>();
    writeVector(_tet,os,&dat);
    writeBinaryData(_weight,os);
    return os.good();
}