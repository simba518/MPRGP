#include "FEMCubatureSolver.h"
#include "FEMUtils.h"
#include "ADMMSolver.h"
#include "solvers/QPSolver.h"
#include <boost/filesystem/operations.hpp>
#include <random>

USE_PRJ_NAMESPACE

//sparse coding solves the following problem:
//argmin_w |Aw-b|_2, |w|_0 <= N s.t. w>=0
//we use matching pursuit greedy solver
template <typename MAT,typename VECB,typename VECX>
scalarD NNLSSmall(const MAT& A,const VECB& b,VECX& x)
{
    Matd ATA=A.transpose()*A;
    Cold ATB=A.transpose()*b;

    x.setZero();
    if(ATA.determinant() == 0.0f)
        return numeric_limits<scalarD>::max();
    else if(A.cols() == 1) {
        x[0]=ATB[0]/ATA(0,0);
        if(x[0] >= 0.0f)
            return (b-A*x).norm();
    } else if(A.cols() == 2) {
        x=ATA.inverse()*ATB;
        if(x[0] >= 0.0f && x[1] >= 0.0f)
            return (b-A*x).norm();
        else {
            x[0]=ATB[0]/ATA(0,0);
            x[1]=ATB[1]/ATA(1,1);
            scalarD res0=numeric_limits<scalarD>::max();
            if(x[0] >= 0.0f)
                res0=(b-A.col(0)*x[0]).norm();
            else x[0]=0.0f;

            scalarD res1=numeric_limits<scalarD>::max();
            if(x[1] >= 0.0f)
                res1=(b-A.col(1)*x[1]).norm();
            else x[1]=0.0f;

            if(res0 < res1) {
                x[1]=0.0f;
                return res0;
            } else if(res1 < res0) {
                x[0]=0.0f;
                return res1;
            }
        }
    } else ASSERT_MSG(false,"Greedy solver only support nr=1/2!")

        x.setZero();
    return b.norm();
}
class SparseCoder
{
public:
    typedef Cold Vec;
    SparseCoder(Matd& A,const Vec& b,sizeType nr)
        :_A(A),_b(b),_nr(nr),_body(NULL),_keyMap(NULL) {
        _scale.resize(_A.cols());
        OMP_PARALLEL_FOR_
        for(sizeType i=0; i<_A.cols(); i++) {
            _scale[i]=1.0f/std::max<scalarD>(_A.col(i).norm(),1E-9f);
            _A.col(i)*=_scale[i];
        }
    }
    void setDebugBody(const FEMBody* body) {
        _body=body;
    }
    void setKeyMap(const vector<sizeType>* keyMap) {
        _keyMap=keyMap;
        ASSERT_MSG(_keyMap->size() == _A.cols()/_nr,"KeyMap size incorrect!");
    }
    const vector<std::pair<Coli,Vec> >& getResult() const {
        return _ret;
    }
    virtual void setParameter(sizeType Ntol,scalarD relTol) =0;
    virtual void solve() =0;
protected:
    scalarD feedResult(const vector<sizeType>& t,const Vec& w,Vec& res) {
        sizeType n=0;
        for(sizeType i=0; i<w.size()/_nr; i++) {
            Eigen::Block<const Vec> BLK=w.block(i*_nr,0,_nr,1);
            if(BLK.norm() > 0.0f)n++;
        }

        std::pair<Coli,Vec> ret;
        ret.first.resize(n);
        ret.second.resize(n*_nr);

        n=0;
        res=_b;
        for(sizeType i=0; i<w.size()/_nr; i++) {
            Eigen::Block<const Vec> BLK=w.block(i*_nr,0,_nr,1);
            Eigen::Block<Vec> scaleBLK=_scale.block(t[i]*_nr,0,_nr,1);
            if(BLK.norm() > 0.0f) {
                ret.first[n]=_keyMap ? (*_keyMap)[t[i]] : t[i];
                ret.second.block(n*_nr,0,_nr,1).array()=BLK.array()*scaleBLK.array();
                res-=_A.block(0,t[i]*_nr,_A.rows(),_nr)*BLK;
                n++;
            }
        }
        if(_body && !_body->_tree.get<bool>("quietCub",false)) {
            std::string path=_body->_tree.get<std::string>("pathCub");
            ostringstream oss;
            oss << path << "/cubature" << ret.first.size() << ".vtk";
            debugWriteVTK(oss.str(),&(ret.first),ret.second);
        }
        _ret.push_back(ret);
        return res.norm()/_b.norm();
    }
    void debugWriteVTK(const std::string& path,Coli* tets,const Vec& weight) const {
        Vec weightAll(_body->nrC());
        weightAll.setZero();
        if(tets) {
            ASSERT_MSG(weight.size() == tets->size()*_nr,"Size mismatch!")
            for(sizeType i=0; i<tets->size(); i++)
                weightAll[(*tets)[i]]=weight.block(i*_nr,0,_nr,1).norm();
        } else {
            ASSERT_MSG(weight.size() == _A.cols(),"Size mismatch!")
            if(_keyMap)ASSERT_MSG(_keyMap->size()*_nr == _A.cols(),"Size mismatch!")
                for(sizeType i=0; i<_A.cols()/_nr; i++) {
                    sizeType t=_keyMap?(*_keyMap)[i]:i;
                    weightAll[t]=weight.block(i*_nr,0,_nr,1).norm();
                }
        }
        _body->writeVTK(VTKWriter<scalar>("CubatureComputation",path,true),NULL,NULL,-1,NULL,&weightAll);
    }
    //data
    Matd& _A;
    Vec _scale;
    const Vec& _b;
    sizeType _nr;
    vector<std::pair<Coli,Vec> > _ret;
    const FEMBody *_body;
    const vector<sizeType>* _keyMap;
};
class GreedySparseCodingSolver : public SparseCoder
{
public:
    GreedySparseCodingSolver(Matd& A,const Vec& b,sizeType nr)
        :SparseCoder(A,b,nr) {
        ASSERT_MSG(nr == 1 || nr == 2,"Greedy solver only support nr=1/2!")
        setParameter(100,1E-2f);
    }
    virtual void setParameter(sizeType Ntol,scalarD tol) {
        _Ntol=Ntol;
        _tol=tol;
        _eps=1E-8f;
        _useMPRGP=false;
    }
    virtual void solve() {
        //output profile data
        std::string path=_body->_tree.get<std::string>("pathCub");
        boost::filesystem::ofstream os(path+"/GreedyCubatureSolverProfile.txt");
        os << "Iter TetId Residual RelRes" << std::endl;

        //temp matrix
        Matd ACurr(_A.rows(),0),ATACurr(0,0);
        Vec ATbCurr;
        Vec wCurr;

        //main loop
        INFO("Greedy Cubature Solver")
        Vec res=_b,weight,L,H;
        sizeType nrC=_A.cols();
        vector<sizeType> ret;
        scalarD relRes;
        weight.setZero(nrC);
        _selected.assign(nrC,false);
        for(sizeType i=0; i<_Ntol; i++) {
            //select next
            ret.push_back(select(res,relRes));

            //setup solver
            sizeType nrT=(sizeType)ret.size();
            if(_useMPRGP) {
                L.setConstant(nrT*_nr,0.0f);
                H.setConstant(nrT*_nr,numeric_limits<scalarD>::max());
                MPRGPQPSolver<scalarD,Kernel<scalarD> > sol;
                sol.reset(L,H);
                sol.setSolverParameters(_eps,10000);
                //sol.setCallback(boost::shared_ptr<Callback<scalarD,Kernel<scalarD> > >(new Callback<scalarD,Kernel<scalarD> >));
                {
                    extend(ACurr,0,_nr);
                    extend(ATACurr,_nr,_nr);
                    extend(ATbCurr,_nr,0);
                    extend(wCurr,_nr,0);
                    for(sizeType t=0,off=(nrT-1)*_nr; t<_nr; t++,off++) {
                        wCurr[off]=0.0f;
                        ACurr.col(off)=_A.col(ret[nrT-1]*_nr+t);
                        ATbCurr[off]=ACurr.col(off).dot(_b);
                        for(sizeType j=0; j<=off; j++)
                            ATACurr(off,j)=ATACurr(j,off)=ACurr.col(off).dot(ACurr.col(j));
                    }
                    //debug code
                    //wCurr.setZero();
                    //ATACurr=ACurr.transpose()*ACurr;
                    //ATbCurr=ACurr.transpose()*_b;
                }
                sol.setMatrix(ATACurr);
                sol.solve(ATbCurr,wCurr);
                ASSERT_MSG(sol.checkKKT(wCurr,ATbCurr),"QP Fail!")
            } else {
                extend(ACurr,0,_nr);
                extend(wCurr,_nr,0);
                for(sizeType t=0,off=(nrT-1)*_nr; t<_nr; t++,off++) {
                    wCurr[off]=0.0f;
                    ACurr.col(off)=_A.col(ret[nrT-1]*_nr+t);
                }
                ASSERT_MSG(NNLS::solve(ACurr,_b,wCurr),"QP Fail!");
            }

            //calculate res and test exit
            scalarD currTol=feedResult(ret,wCurr.array(),res);
            INFOV("tets: %d, tol: %f, relRes: %f",(sizeType)ret.size(),currTol,relRes)
            os << (i+1) << " " << ret[ret.size()-1] << " " << currTol << " " << relRes << std::endl;
            if(currTol < _tol)
                break;
        }
    }
    sizeType select(const Vec& res,scalarD& relRes) const {
        //we solve a sub-problem:
        //min|A*w-res|_2
        //s.t.w>=0
        Vec x(_nr);
        sizeType minId=-1;
        scalarD minVal=numeric_limits<scalarD>::max();
        for(sizeType i=0; i<_A.cols()/_nr; i++)
            if(!_selected[i]) {
                scalarD resI=NNLSSmall(_A.block(0,i*_nr,_A.rows(),_nr),res,x);
                if(resI < minVal) {
                    minId=i;
                    minVal=resI;
                }
            }
        relRes=minVal/res.norm();
        return minId;
    }
protected:
    //data
    sizeType _Ntol;
    scalarD _tol,_eps;
    //temporary data
    vector<bool> _selected;
    bool _useMPRGP;
};
class ReweightedL0SparseCodingSolver : public SparseCoder
{
public:
    ReweightedL0SparseCodingSolver(Matd& A,const Vec& b,sizeType nr)
        :SparseCoder(A,b,nr) {
        _sol.reset(new ADMMSolver<Kernel<scalarD> >(_A));
    }
    virtual void setParameter(sizeType Ntol,scalarD relTol) {
        _relTol=relTol;
        _sqrDelta=_b.squaredNorm()*_relTol*_relTol;
        _p=0.3f;
        _thres=1E-5f;
    }
    virtual void solve() {
        //initialize
        scalarD currP=1.0f;
        sizeType nrX=_A.cols();
        //output profile data
        std::string path=_body->_tree.get<std::string>("pathCub");
        boost::filesystem::ofstream os(path+"/L0CubatureSolverProfile.txt");

        //ADMM for L_{_p}/L2con+ optimization
        Vec x=Vec::Constant(nrX,_thres*2.0f);
        for(sizeType i=0; currP > _p; i++) {
            _sol->setProb(x,_b,currP,_thres);
            ASSERT_MSG(_sol->solve(_sqrDelta),"L1 Fail!");
            _sol->checkBound(_sqrDelta);
            _sol->fetchX(x);
            if(_body && !_body->_tree.get<bool>("quietCub",false)) {
                ostringstream oss;
                oss << path << "/L1weight" << i << ".vtk";
                debugWriteVTK(oss.str(),NULL,x);
            }
            currP-=0.1f;
        }

        //round result x
        scalarD clampThres=0.01f*x.maxCoeff();
        vector<sizeType> ret;
        for(sizeType i=0; i<nrX; i+=_nr)
            if(x.block(i,0,_nr,1).cwiseAbs().maxCoeff() > clampThres)
                ret.push_back(i/_nr);

        //run QP to find weight
        sizeType nrT=(sizeType)ret.size();
        Matd ABlk(_A.rows(),nrT*_nr);
        for(sizeType i=0; i<nrT; i++)
            ABlk.block(0,i*_nr,_A.rows(),_nr)=
                _A.block(0,ret[i]*_nr,_A.rows(),_nr);
        Vec weight(nrT*_nr);
        NNLS::solve(ABlk,_b,weight);

        Vec res;
        scalarD relTol=feedResult(ret,weight,res);
        INFOV("Nr Tet: %d, RelTol: %f",(sizeType)ret.size(),relTol)
        os << (sizeType)ret.size() << " " << relTol;
    }
protected:
    //data
    boost::shared_ptr<ADMMSolver<Kernel<scalarD> > > _sol;
    scalarD _relTol,_sqrDelta,_p,_thres;
};

//FEM cubature solver
void FEMCubatureProb::transformMetric(const Eigen::SparseMatrix<scalarD,0,sizeType>& weight,const Matd& B,Matd& M)
{
    Matd B1=B;
    EigenSolver::makeOrthogonal(&weight,B1);

    Matd invB1TB1=(B1.transpose()*B1).eval().inverse();
    M=invB1TB1*(B1.transpose()*B)*M;
}
FEMCubatureSolver::FEMCubatureSolver(FEMBody& body,FEMCubatureProb& obj)
    :_body(body),_obj(obj)
{
    if(_body._tree.find("pathCub") == _body._tree.not_found())
        _body._tree.put<std::string>("pathCub","meshCubature");
    if(_body._tree.find("overwriteCub") == _body._tree.not_found())
        _body._tree.put<bool>("overwriteCub",false);
    sizeType N=std::max<sizeType>(std::min<sizeType>(body.nrC()/10,100),1);
    if(_body._tree.find("nrCub") == _body._tree.not_found())
        _body._tree.put<sizeType>("nrCub",N);
    if(_body._tree.find("nrPose") == _body._tree.not_found())
        _body._tree.put<sizeType>("nrPose",std::min<sizeType>(N*10,1000));
    if(_body._tree.find("tolSC") == _body._tree.not_found())
        _body._tree.put<scalar>("tolSC",0.01f);
    INFOV("Cubature Param: %d cells %d poses %f relErr",
          _body._tree.get<sizeType>("nrCub"),
          _body._tree.get<sizeType>("nrPose"),
          _body._tree.get<scalar>("tolSC"))
}
void FEMCubatureSolver::generate(const Vec& weight,bool debugWrite)
{
    //identify LMA basis
    //we use random sampling for training data
    //we have to use LMA for this sampling process
    const Matd& BASIS=_body._basis->_U;
    sizeType nrLMA=_body._basis->nrLMA(BASIS);
    Vec UTKUD=(BASIS.transpose()*(_body._basis->_K*BASIS).eval()).diagonal();
    Vec UTMUD=(BASIS.transpose()*(_body._basis->_M*BASIS).eval()).diagonal();

    //find approximately the coefficient for sample scale
    scalarD coef;
    {
        BBox<scalar> bb;
        for(sizeType i=0; i<_body.nrV(); i++)
            bb.setUnion(_body.getV(i)._pos0);
        scalar dist=bb.getExtent().norm();
        Vec b0=BASIS.col(0).cwiseAbs();
        coef=(dist*0.1f/(b0.array()+1E-6f)).matrix().minCoeff()*UTKUD[0];
    }

    //generate training data
    std::string path=_body._tree.get<std::string>("pathCub");
    bool overwrite=_body._tree.get<bool>("overwriteCub");
    sizeType nrPose=_body._tree.get<sizeType>("nrPose");
    sizeType nrC=_body.nrC();
    sizeType rows=BASIS.rows();
    Vec2i nr=_obj.nr();
    if(boost::filesystem::exists("./"+path+"/TrainingData.dat")) {
        if(!read(boost::filesystem::ifstream("./"+path+"/TrainingData.dat",ios::binary)))
            overwrite=true;
        else
            overwrite=overwrite ||
                      _szA[0]%nr[0] != 0 || _szA[1] != nrC*nr[1] ||
                      _b.size() != _szA[0] || _coef.size() != _b.size() ||
                      _S.rows() != nrLMA || _S.cols() != nrPose;
    } else overwrite=true;

    //size incorrect or insufficient training data
    if(overwrite) {
        boost::filesystem::create_directory(path);
        if(debugWrite)
            boost::filesystem::create_directory("./"+path+"/TrainingPose");

        _szA[0]=nrPose*nr[0];
        _szA[1]=nrC*nr[1];
        Matd ABlk(_szA[0],nr[1]);
        _b.setZero(nrPose*nr[0]);
        _w=weight;
        _coef.setZero(nrPose*nr[0]);
        _S.setZero(nrLMA,nrPose);
        INFO("Generating Pose!")
        Matd U=BASIS.block(0,0,rows,nrLMA);
        std::default_random_engine generator;
        for(sizeType i=0; i<nrPose; i++) {
            for(sizeType s=0; s<nrLMA; s++)
                _S(s,i)=std::normal_distribution<scalarD>(0.0f,coef/UTKUD[s])(generator);
            if(debugWrite) {
                _body.setDPos(U*_S.col(i));
                ostringstream oss;
                oss << "./"+path+"/TrainingPose/pose" << i << ".vtk";
                VTKWriter<scalar> f("Pose",oss.str(),true);
                _body.writeVTK(f);
            }
        }

        boost::filesystem::fstream os("./"+path+"/TrainingData.dat",ios::binary|ios::in|ios::out|ios::trunc);
        writeBinaryData(_szA[0],os);
        writeBinaryData(_szA[1],os);
        INFO("Generating Prob!")
        ASSERT_MSGV(weight.size() == _szA[1],"Size of weight mismatch!")
        //generate _A
        scalarD Anorm=0.0f;
        for(sizeType e=0; e<nrC; e++) {
            if(e > 0 && e%1000 == 0)
                INFOV("Generated %d cells!",e)
                ABlk=_obj(e,_S);
            Anorm+=ABlk.squaredNorm();
            os.write((char*)ABlk.data(),ABlk.size()*sizeof(Matd::Scalar));
            _b+=ABlk*_w.block(e*nr[1],0,nr[1],1);
        }
        INFOV("ANorm: %f, bNorm: %f",sqrt(Anorm),_b.norm())
        _obj.debugCallback(_b,_S);
        //normalize _b at each pose
        OMP_PARALLEL_FOR_
        for(sizeType i=0; i<nrPose; i++) {
            scalarD coef=1.0f/_b.block(i*nr[0],0,nr[0],1).norm();
            _coef.block(i*nr[0],0,nr[0],1).setConstant(coef);
            _b.block(i*nr[0],0,nr[0],1)*=coef;
        }
        writeBinaryData(_b,os);
        writeBinaryData(_w,os);
        writeBinaryData(_coef,os);
        writeBinaryData(_S,os);
        scaleCol(os);
    }
}
void FEMCubatureSolver::solve(boost::shared_ptr<Cubature>& ptrCub,const std::string& extension)
{
    if(!ptrCub)
        ptrCub.reset(new Cubature);
    Cubature& cub=*ptrCub;

    //test existing result
    std::string path=_body._tree.get<std::string>("pathCub");
    if(!_body._tree.get<bool>("overwriteCub")) {
        ostringstream oss;
        oss << path << "/mesh" << extension;

        Cubature cubFile;
        boost::filesystem::ifstream f(oss.str(),ios::binary);
        if(boost::filesystem::exists(oss.str()) && cubFile.read(f)) {
            cub=cubFile;
            return;
        }
    }

    //if not, solve sparse coding
    _body.clearPos();
    vector<std::pair<Coli,Vec> > rets;
    solveSparseCodingHier(rets);
    boost::filesystem::create_directory(path);
    for(sizeType i=0; i<(sizeType)rets.size(); i++) {
        sizeType nrT=(sizeType)rets[i].first.size();
        sizeType dim=_body.dim();

        cub._weight=rets[i].second;
        cub._tet.resize(nrT*(dim+1));
        for(sizeType j=0; j<nrT; j++)
            for(char v=0; v<(dim+1); v++) {
                cub._tet[j*(dim+1)+v]._id=rets[i].first[j];
                cub._tet[j*(dim+1)+v]._coef=Vec4::Unit(v);
            }
        {
            ostringstream oss;
            oss << path << "/mesh" << nrT << extension;
            boost::filesystem::ofstream of(oss.str(),ios::binary);
            cub.write(of);
        }
    }
    {
        ostringstream oss;
        oss << path << "/mesh" << extension;
        boost::filesystem::ofstream of(oss.str(),ios::binary);
        cub.write(of);
    }
}
void FEMCubatureSolver::solveSparseCodingHier(vector<std::pair<Coli,Vec> >& ret)
{
    bool quiet=_body._tree.get<bool>("quietCub",false);
    std::string path=_body._tree.get<std::string>("pathCub");
    sizeType nrP=_body._tree.get<sizeType>("hierThres",3000);
    boost::filesystem::ifstream is("./"+path+"/TrainingData.dat",ios::binary);
    read(is);
    if(_body.nrC() > nrP) {
        vector<sizeType> keyMap;
        //we observed that approximately
        //nrP*3 particles are returned
        nrP/=3;
        SegmentBody::segment(_body,nrP,keyMap);
        if(!quiet)
            SegmentBody::writeVTK(path+"/group.vtk",_body,keyMap);

        //form coarse sparse coding
        Matd AGroup=Matd::Zero(_szA[0],nrP);
        groupCol(AGroup,is,keyMap,nrP);

        //solve coarse sparse coding
        _body._tree.put<bool>("quietCub",true);
        solveSparseCoding(AGroup,ret,1,NULL);
        Coli groups=ret.back().first;
        boost::unordered_set<sizeType> groupSet;
        for(sizeType i=0; i<groups.size(); i++)
            groupSet.insert(groups[i]);

        //find invKeyMap
        vector<sizeType> invKeyMap;
        for(sizeType i=0; i<_body.nrC(); i++)
            if(groupSet.find(keyMap[i]) != groupSet.end())
                invKeyMap.push_back(i);
        if(!quiet)
            SegmentBody::writeSelVTK(path+"/selGroup.vtk",_body,invKeyMap);

        //form fine sparse coding
        selectCol(AGroup,is,&invKeyMap);
        _body._tree.put<bool>("quietCub",quiet);
        solveSparseCoding(AGroup,ret,_obj.nr()[1],&invKeyMap);
    } else {
        _body._tree.put<bool>("quietCub",quiet);
        Matd A;
        selectCol(A,is,NULL);
        solveSparseCoding(A,ret,_obj.nr()[1],NULL);
    }
}
void FEMCubatureSolver::solveSparseCoding(Matd& A,vector<std::pair<Coli,Vec> >& ret,sizeType nr,vector<sizeType>* keyMap) const
{
    boost::shared_ptr<SparseCoder> sol;
    sizeType option=_body._tree.get<sizeType>("optionCub",2);
    if(option==1)
        sol.reset(new GreedySparseCodingSolver(A,_b,nr));
    else if(option==2)
        sol.reset(new ReweightedL0SparseCodingSolver(A,_b,nr));
    else ASSERT_MSG(false,"Unknown Cubature Solver Type!")
        if(keyMap)sol->setKeyMap(keyMap);
    sol->setDebugBody(&_body);

    sizeType nrCub=_body._tree.get<sizeType>("nrCub");
    scalarD tolSC=_body._tree.get<scalarD>("tolSC");
    sol->setParameter(nrCub,tolSC);
    sol->solve();
    ret=sol->getResult();
}
const Matd& FEMCubatureSolver::getS() const
{
    return _S;
}
//handle out-of-core A
bool FEMCubatureSolver::read(std::istream& is)
{
    if(!readBinaryData(_szA[0],is).good())return false;
    if(!readBinaryData(_szA[1],is).good())return false;

    is.seekg(_szA[0]*_szA[1]*sizeof(Matd::Scalar)+sizeof(sizeType)*2,ios::beg);
    if(!is.good())return false;

    if(!readBinaryData(_b,is).good())return false;
    if(!readBinaryData(_w,is).good())return false;
    if(!readBinaryData(_coef,is).good())return false;
    if(!readBinaryData(_S,is).good())return false;

    sizeType nrPose=_body._tree.get<sizeType>("nrPose");
    if(_szA[0] != _b.size() || _szA[1] != _w.size() ||
            _szA[0] != _coef.size() || nrPose != _S.cols())
        return false;
    return true;
}
void FEMCubatureSolver::scaleCol(std::iostream& ss) const
{
    sizeType off=sizeof(sizeType)*2;
    sizeType szC=sizeof(Matd::Scalar)*_szA[0];

    Vec ABlk(_szA[0]);
    for(sizeType i=0; i<_szA[1]; i++) {
        ss.seekg(off+szC*i,ios_base::beg);
        ss.read((char*)ABlk.data(),szC);
        ABlk.array()*=_coef.array();

        ss.seekp(off+szC*i,ios_base::beg);
        ss.write((char*)ABlk.data(),szC);
    }
}
void FEMCubatureSolver::groupCol(Matd& AGroup,std::istream& is,const vector<sizeType>& keyMap,sizeType nrG) const
{
    sizeType off=sizeof(sizeType)*2;
    sizeType szC=sizeof(Matd::Scalar)*_szA[0];

    sizeType nr=_obj.nr()[1];
    sizeType nrC=_body.nrC();
    Matd cache(_szA[0],nr);
    AGroup.setZero(_szA[0],nrG);
    for(sizeType i=0; i<nrC; i++) {
        is.seekg(off+szC*i*nr,ios_base::beg);
        is.read((char*)cache.data(),szC*nr);
        AGroup.col(keyMap[i])+=cache*_w.block(i*nr,0,nr,1);
    }
}
void FEMCubatureSolver::selectCol(Matd& AGroup,std::istream& is,const vector<sizeType>* keyMap) const
{
    sizeType off=sizeof(sizeType)*2;
    sizeType szC=sizeof(Matd::Scalar)*_szA[0];

    sizeType nr=_obj.nr()[1];
    sizeType nrC=keyMap?keyMap->size():_body.nrC();
    AGroup.resize(_szA[0],nrC*nr);
    for(sizeType i=0; i<nrC; i++) {
        sizeType c=keyMap?(*keyMap)[i]:i;
        is.seekg(off+szC*c*nr,ios_base::beg);
        is.read((char*)(AGroup.data()+_szA[0]*i*nr),szC*nr);
    }
}

//CFEMBody
CFEMBody::CFEMBody() {}
CFEMBody::CFEMBody(const FEMBody& body,const Cubature& ce)
{
    //build a body representing all cubatures
    sizeType dim=body.dim();
    sizeType nrE=(sizeType)ce._tet.size()/(dim+1);

    _vss.resize(nrE*(dim+1));
    _css.resize(nrE);
    for(sizeType i=0; i<nrE; i++) {
        _css[i].reset(new FEMCell);
        FEMCell& c=*(_css[i]);
        for(char v=0; v<dim+1; v++) {
            sizeType off=i*(dim+1)+v;
            _vss[off].reset(new FEMVertex);
            c._v[v]=_vss[off];
        }
    }

    //build interpolation matrix
    sizeType nrV=(sizeType)ce._tet.size();
    vector<Eigen::Triplet<scalarD,sizeType> > H;
    for(sizeType i=0; i<nrV; i++) {
        const FEMInterp& I=ce._tet[i];
        for(sizeType vi=0; vi<body.dim()+1; vi++)
            addI3x3(H,i*3,body.getC(I._id)._v[vi]->_index*3,I._coef[vi]);
    }
    _cInterp.resize(nrV*3,body.nrV()*3);
    _cInterp.setFromTriplets(H.begin(),H.end());

    //assemble
    Vec pos;
    body.getPos(pos);

    Vec pos0=_cInterp*pos;
    for(sizeType i=0; i<nrV; i++)
        getV(i)._pos0=pos0.block<3,1>(i*3,0).cast<scalar>();
    assemble();
}
FEMBody& CFEMBody::operator=(const FEMBody& other)
{
    const CFEMBody* otherDef=static_cast<const CFEMBody*>(&other);
    ASSERT_MSG(otherDef,"Different Body Type!");
    FEMBody::operator=(other);
    _cInterp=otherDef->_cInterp;
    return *this;
}
boost::shared_ptr<Serializable> CFEMBody::copy() const
{
    return boost::shared_ptr<Serializable>(new CFEMBody);
}