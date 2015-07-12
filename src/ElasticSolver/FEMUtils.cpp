#include "FEMUtils.h"
#include "FEMMesh.h"
#include "FEMBVHBuilder.h"
#include "FEMSparseReducedBasis.h"
#include "ParallelPoissonDiskSampling.h"
#include "solvers/AINVPreconditioner.h"
#include "solvers/DACGSolver.h"
#include <boost/unordered_set.hpp>
extern "C"
{
    /* Subroutine */
    int nnls_(double *a,int *mda,int *m,int *n,
              double *b,double *x,double *rnorm,double *w,double *zz,
              int *index,int *mode);
}

USE_PRJ_NAMESPACE

size_t Hash::operator()(const Vec2i& key) const
{
    boost::hash<sizeType> h;
    return h(key[0])+h(key[1]);
}
size_t Hash::operator()(const Vec3i& key) const
{
    boost::hash<sizeType> h;
    return h(key[0])+h(key[1])+h(key[2]);
}
size_t Hash::operator()(const Vec4i& key) const
{
    boost::hash<sizeType> h;
    return h(key[0])+h(key[1])+h(key[2])+h(key[3]);
}
//mesh segmentation
struct KeyPointCallback {
    KeyPointCallback(sizeType& id):_id(id) {}
    void updateDist(const Node<sizeType>& node,const Vec3& pt,Vec3& cb,Vec3& n,scalar& dist) const {
        scalar newDist=(pt-node._bb._minC).norm();
        if(node._nrCell==1 && newDist < dist) {
            dist=newDist;
            _id=node._cell;
        }
    }
    scalar depth() const {
        return numeric_limits<scalar>::max();
    }
    sizeType& _id;
};
void SegmentBody::segment(const FEMBody& body,sizeType& nrP,vector<sizeType>& keyMap,bool debugPt)
{
    //use two level sparse coding
    vector<Node<sizeType> > keyBVH;
    {
        //sample a set of points
        ParallelPoissonDiskSampling sampler(body.dim());
        //vss
        sizeType nrV=body.nrV();
        std::vector<Vec3,Eigen::aligned_allocator<Vec3> > vss(nrV);
        OMP_PARALLEL_FOR_
        for(sizeType i=0; i<nrV; i++)
            vss[i]=body.getV(i)._pos0;
        //css
        sizeType nrC=body.nrC();
        std::vector<Vec4i,Eigen::aligned_allocator<Vec4i> > css(nrC,Vec4i::Constant(-1));
        OMP_PARALLEL_FOR_
        for(sizeType i=0; i<nrC; i++)
            for(char v=0; v<(body.dim()+1); v++)
                css[i][v]=body.getC(i)._v[v]->_index;
        sampler.sample(vss,css,&nrP);
        if(debugPt)
            sampler.getPSet().writeVTK("./samplePSet.vtk");
        //insert all keys
        nrP=sampler.getPSet().size();
        keyBVH.resize(nrP);
        OMP_PARALLEL_FOR_
        for(sizeType i=0; i<nrP; i++) {
            const ParticleN<scalar>& p=sampler.getPSet()[i];
            keyBVH[i]._nrCell=1;
            keyBVH[i]._bb=BBox<scalar>(p._pos,p._pos);
            keyBVH[i]._cell=i;
        }
        //build a BVH for this particle set
        if(body.dim() == 2) {
            BVHBuilder<Node<sizeType>,2> builder;
            builder.buildBVH(keyBVH);
        } else {
            BVHBuilder<Node<sizeType>,3> builder;
            builder.buildBVH(keyBVH);
        }
        OMP_PARALLEL_FOR_
        for(sizeType i=nrP; i<(sizeType)keyBVH.size(); i++)
            keyBVH[i]._cell=-1;
    }

    //find the map
    keyMap.resize(body.nrC());
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<body.nrC(); i++) {
        scalar dist=numeric_limits<scalar>::max();
        Vec3 cp,n,ctr=Vec3::Zero();
        for(char v=0; v<(body.dim()+1); v++)
            ctr+=body.getC(i)._v[v]->_pos0;
        ctr/=(scalar)(body.dim()+1);
        KeyPointCallback cb(keyMap[i]);
        BVHQuery<sizeType>(keyBVH,body.dim(),-1).pointDistQuery(ctr,cb,cp,n,dist);
    }
}
void SegmentBody::writeVTK(const std::string& path,const FEMBody& body,const vector<sizeType>& keyMap)
{
    Cold css(body.nrC());
    for(sizeType i=0; i<body.nrC(); i++)
        css[i]=(scalarD)keyMap[i];
	VTKWriter<scalar> w("SegmentResult",path,true);
    body.writeVTK(w,NULL,NULL,-1,NULL,&css);
}
void SegmentBody::writeSelVTK(const std::string& path,const FEMBody& body,const vector<sizeType>& keyMap)
{
    Cold css=Cold::Zero(body.nrC());
    for(sizeType i=0; i<(sizeType)keyMap.size(); i++)
        css[keyMap[i]]=1.0f;
	VTKWriter<scalar> w("SegmentResult",path,true);
    body.writeVTK(w,NULL,NULL,-1,NULL,&css);
}
//nonnegative least square
bool NNLS::solve(const Matd& A,const Cold& b,Cold& x)
{
    Matd Atmp=A;
    int mda=(int)Atmp.rows();
    int m=(int)Atmp.rows();
    int n=(int)Atmp.cols();
    Cold btmp=b;
    double rnorm;
    Cold w(n);
    Cold zz(m);
    Eigen::Matrix<int,-1,1> index(n);
    int mode;
    nnls_(Atmp.data(),&mda,&m,&n,
          btmp.data(),x.data(),&rnorm,w.data(),zz.data(),
          index.data(),&mode);
    return mode == 1;
}
//basis related
class ModeBasisSolver
{
public:
    ModeBasisSolver(const Eigen::SparseMatrix<scalarD,0,sizeType>& K,const Eigen::SparseMatrix<scalarD,0,sizeType>& M,const boost::unordered_set<sizeType>& fixSet) {
        int j=0;
        _indexMap.assign(K.rows(),-1);
        for(int i=0; i<K.rows(); i++)
            if(fixSet.find(i) == fixSet.end())
                _indexMap[i]=j++;
        _compressedSize=j;
        _factorize=true;
        compress(K,_K);
        compress(M,_M);
    }
    void solve(const Eigen::Matrix<scalarD,-1,1>& pos,Eigen::Matrix<scalarD,-1,-1>& U,Eigen::Matrix<scalarD,-1,1>& lambda,sizeType excludeRigid) {
        Eigen::Matrix<scalarD,-1,-1> UTmp;
        solveDACG(pos,UTmp,lambda,excludeRigid);
        extend(UTmp,U);
    }
    void setFactorize(bool factorize) {
        _factorize=factorize;
    }
protected:
    void solveDACG(const Eigen::Matrix<scalarD,-1,1>& pos,Eigen::Matrix<scalarD,-1,-1>& U,Eigen::Matrix<scalarD,-1,1>& lambda,sizeType excludeRigid) {
        sizeType n=(sizeType)lambda.size();
        lambda.setRandom();
        U.resize(_K.rows(),n);
        U.setRandom();

        Matd U0(_K.rows(),0);
        if(excludeRigid > 0) {
            if(excludeRigid == 2)
                U0=build2DRigidMode(pos);
            else U0=build3DRigidMode(pos);
            INFOV("Kernel Space: %f",(_K*U0).norm());
        }

        if(_factorize) {
            DACGSolver<scalarD,Kernel<scalarD>,DirectPreconSolver<scalarD> > solver;
            solver.setCallback(boost::shared_ptr<Callback<scalarD,Kernel<scalarD> > >(new Callback<scalarD,Kernel<scalarD> >));
            solver.setSolverParameters(1E-20f,1E-5f,1000000);
            solver.setKrylovA(_K);
            solver.setKrylovB(_M);
            //set preconditioner
            _KPre=_K;
            EigenSolver::scaleDiagonal(_KPre,1.01f,1E-2f);
            static_cast<DirectPreconSolver<scalarD>*>(solver.getPre())->setMatrix(_KPre,true);
            solver.setU0(U0);
            solver.solve(n,lambda,U);
        } else {
            DACGSolver<scalarD,Kernel<scalarD>,SymAINVPreconSolver<scalarD> > solver;
            solver.setCallback(boost::shared_ptr<Callback<scalarD,Kernel<scalarD> > >(new Callback<scalarD,Kernel<scalarD> >));
            solver.setSolverParameters(1E-20f,1E-5f,1000000);
            solver.setKrylovA(_K);
            solver.setKrylovB(_M);
            //set preconditioner
            _KPre=_K;
            EigenSolver::scaleDiagonal(_KPre,1.01f,1E-2f);
            static_cast<SymAINVPreconSolver<scalarD>*>(solver.getPre())->setSolverParameter(0.01f,0.01f);
            static_cast<SymAINVPreconSolver<scalarD>*>(solver.getPre())->setMatrix(_KPre,true);
            solver.setU0(U0);
            solver.solve(n,lambda,U);
        }
    }
    Matd build2DRigidMode(const Eigen::Matrix<scalarD,-1,1>& pos) const {
        Matd ret(_K.rows(),3);
        ret.setZero();
        //translation
        for(int i=0; i<ret.rows(); i+=2) {
            ret.block<2,2>(i,0).setIdentity();
            Vec2d posI=pos.block<2,1>(i*3/2,0);
            ret.block<2,1>(i,2)=Vec2d(posI[1],-posI[0]);
        }
        return ret;
    }
    Matd build3DRigidMode(const Eigen::Matrix<scalarD,-1,1>& pos) const {
        Matd ret(_K.rows(),6);
        ret.setZero();
        //translation
        for(int i=0; i<ret.rows(); i+=3) {
            ret.block<3,3>(i,0).setIdentity();
            Vec3d posI=pos.block<3,1>(i,0);
            ret.block<3,1>(i,3)=Vec3d(    0.0f, posI[2],-posI[1]);
            ret.block<3,1>(i,4)=Vec3d(-posI[2],    0.0f, posI[0]);
            ret.block<3,1>(i,5)=Vec3d( posI[1],-posI[0],    0.0f);
        }
        return ret;
    }
    void compress(const Eigen::SparseMatrix<scalarD,0,sizeType>& from,Eigen::SparseMatrix<scalarD,0,sizeType>& to) const {
        vector<Eigen::Triplet<scalarD,sizeType> > H;
        for(int k=0; k<from.outerSize(); ++k)
            for(Eigen::SparseMatrix<scalarD,0,sizeType>::InnerIterator it(from,k); it; ++it)
                if(_indexMap[it.row()] != -1 && _indexMap[it.col()] != -1)
                    H.push_back(Eigen::Triplet<scalarD,sizeType>(_indexMap[it.row()],_indexMap[it.col()],it.value()));
        to.resize(_compressedSize,_compressedSize);
        to.setFromTriplets(H.begin(),H.end());
    }
    void extend(const Eigen::Matrix<scalarD,-1,-1>& from,Eigen::Matrix<scalarD,-1,-1>& to) const {
        to.setZero(_indexMap.size(),from.cols());
        for(sizeType i=0; i<(sizeType)_indexMap.size(); i++)
            if(_indexMap[i] != -1)
                to.row(i)=from.row(_indexMap[i]);
    }
    //data
    bool _factorize;
    sizeType _compressedSize;
    vector<sizeType> _indexMap;
    Eigen::SparseMatrix<scalarD,0,sizeType> _K,_KPre,_M;
};
void EigenSolver::solveEigen(const FEMBody& body,const boost::unordered_set<sizeType>& fixSet,Eigen::Matrix<scalarD,-1,1>& lambda,bool excludeRigid)
{
    SparseReducedBasis& basis=*(body._basis);
    ModeBasisSolver sol(basis._K,basis._M,fixSet);
    sol.setFactorize(body._tree.get<bool>("useDirectFactorize",true));

    Cold pos;
    body.getPos(pos);
    sol.solve(pos,basis._U,lambda,excludeRigid?body.dim():0);
}
void EigenSolver::makeOrthogonal(const Eigen::SparseMatrix<scalarD,0,sizeType>* A,Matd& U)
{
#define MUL(II) (A?((*A)*II).eval():II)
    sizeType n=U.cols();
    for(int i=0; i<n; i++) {
        for(int j=0; j<i; j++) {
            const scalarD alpha=U.col(i).dot(MUL(U.col(j)))/
                                U.col(j).dot(MUL(U.col(j)));
            U.col(i)-=alpha*U.col(j);
        }
        U.col(i)/=sqrt(U.col(i).dot(MUL(U.col(i))));
    }
#undef MUL
}
void EigenSolver::scaleDiagonal(Eigen::SparseMatrix<scalarD,0,sizeType>& A,scalarD coef,scalarD minV)
{
    //add diagonal regularizer
    for(sizeType i=0; i<A.rows(); i++) {
        scalarD& diagV=A.coeffRef(i,i);
        diagV=std::max(diagV*coef,minV);
    }
}
