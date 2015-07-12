#include "ClothCollision.h"
// #include <optimization.h>
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

ClothCollision::CollisionHandler::CollisionHandler():_id(0) {}
void ClothCollision::CollisionHandler::handle(boost::shared_ptr<ClothMesh::ClothVertex> V1,boost::shared_ptr<ClothMesh::ClothTriangle> T2,const Vec3d n,const Vec4d& omg,scalarD t)
{
    const ClothMesh::ClothVertex& v0=*V1;
    const ClothMesh::ClothVertex& v1=*(T2->getV0());
    const ClothMesh::ClothVertex& v2=*(T2->getV1());
    const ClothMesh::ClothVertex& v3=*(T2->getV2());

    vector<scalarD> css;
    vector<Vec3d,Eigen::aligned_allocator<Vec3d> > tvV;
    tvV.push_back(v0._lastPos);
    tvV.push_back(v1._lastPos);
    tvV.push_back(v2._lastPos);
    tvV.push_back(v3._lastPos);
    tvV.push_back(interp1D(v0._lastPos,v0._pos,t));
    tvV.push_back(interp1D(v1._lastPos,v1._pos,t));
    tvV.push_back(interp1D(v2._lastPos,v2._pos,t));
    tvV.push_back(interp1D(v3._lastPos,v3._pos,t));
    tvV.push_back(tvV[4]+n);
    tvV.push_back(tvV[5]-n);
    tvV.push_back(tvV[6]-n);
    tvV.push_back(tvV[7]-n);

    boost::filesystem::create_directory("./coll");
    ostringstream oss;
    oss << "./coll/frm" << _id << ".vtk";
    VTKWriter<scalarD> os("coll",oss.str(),true);
    os.appendPoints(tvV.begin(),tvV.end());

    vector<Vec3i,Eigen::aligned_allocator<Vec3i> > tvI;
    tvI.push_back(Vec3i::Constant(0));
    tvI.push_back(Vec3i::Constant(4));
    css.push_back(0.0f);
    css.push_back(0.0f);
    os.appendCells(tvI.begin(),tvI.end(),VTKWriter<scalarD>::POINT);

    tvI.clear();
    tvI.push_back(Vec3i(1,2,3));
    tvI.push_back(Vec3i(5,6,7));
    css.push_back(1.0f);
    css.push_back(1.0f);
    os.appendCells(tvI.begin(),tvI.end(),VTKWriter<scalarD>::TRIANGLE);

    tvI.clear();
    tvI.push_back(Vec3i::Constant(0));
    tvI.push_back(Vec3i::Constant(1));
    tvI.push_back(Vec3i::Constant(2));
    tvI.push_back(Vec3i::Constant(3));
    tvI.push_back(Vec3i::Constant(4));
    tvI.push_back(Vec3i::Constant(5));
    tvI.push_back(Vec3i::Constant(6));
    tvI.push_back(Vec3i::Constant(7));
    css.push_back(2.0f);
    css.push_back(2.0f);
    css.push_back(2.0f);
    css.push_back(2.0f);
    css.push_back(2.0f);
    css.push_back(2.0f);
    css.push_back(2.0f);
    css.push_back(2.0f);
    os.appendCells(tvI.begin(),tvI.end(),VTKWriter<scalarD>::POINT);

    tvI.clear();
    tvI.push_back(Vec3i(4,8,0));
    tvI.push_back(Vec3i(5,9,0));
    tvI.push_back(Vec3i(6,10,0));
    tvI.push_back(Vec3i(7,11,0));
    css.push_back(3.0f);
    css.push_back(3.0f);
    css.push_back(3.0f);
    css.push_back(3.0f);
    os.appendCells(tvI.begin(),tvI.end(),VTKWriter<scalarD>::LINE);
    os.appendCustomData("Color",css.begin(),css.end());
    _id++;
}
void ClothCollision::CollisionHandler::handle(boost::shared_ptr<ClothMesh::ClothEdge> E1,boost::shared_ptr<ClothMesh::ClothEdge> E2,const Vec3d n,const Vec4d& omg,scalarD t)
{
    const ClothMesh::ClothVertex& v0=*(E1->_v[0]);
    const ClothMesh::ClothVertex& v1=*(E1->_v[1]);
    const ClothMesh::ClothVertex& v2=*(E2->_v[0]);
    const ClothMesh::ClothVertex& v3=*(E2->_v[1]);

    vector<scalarD> css;
    vector<Vec3d,Eigen::aligned_allocator<Vec3d> > tvV;
    tvV.push_back(v0._lastPos);
    tvV.push_back(v1._lastPos);
    tvV.push_back(v2._lastPos);
    tvV.push_back(v3._lastPos);
    tvV.push_back(interp1D(v0._lastPos,v0._pos,t));
    tvV.push_back(interp1D(v1._lastPos,v1._pos,t));
    tvV.push_back(interp1D(v2._lastPos,v2._pos,t));
    tvV.push_back(interp1D(v3._lastPos,v3._pos,t));
    tvV.push_back(tvV[4]+n);
    tvV.push_back(tvV[5]+n);
    tvV.push_back(tvV[6]-n);
    tvV.push_back(tvV[7]-n);

    boost::filesystem::create_directory("./coll");
    ostringstream oss;
    oss << "./coll/frm" << _id << ".vtk";
    VTKWriter<scalarD> os("coll",oss.str(),true);
    os.appendPoints(tvV.begin(),tvV.end());

    vector<Vec3i,Eigen::aligned_allocator<Vec3i> > tvI;
    tvI.push_back(Vec3i(0,1,0));
    tvI.push_back(Vec3i(2,3,0));
    tvI.push_back(Vec3i(4,5,0));
    tvI.push_back(Vec3i(6,7,0));
    css.push_back(0.0f);
    css.push_back(1.0f);
    css.push_back(0.0f);
    css.push_back(1.0f);
    os.appendCells(tvI.begin(),tvI.end(),VTKWriter<scalarD>::LINE);

    tvI.clear();
    tvI.push_back(Vec3i::Constant(0));
    tvI.push_back(Vec3i::Constant(1));
    tvI.push_back(Vec3i::Constant(2));
    tvI.push_back(Vec3i::Constant(3));
    tvI.push_back(Vec3i::Constant(4));
    tvI.push_back(Vec3i::Constant(5));
    tvI.push_back(Vec3i::Constant(6));
    tvI.push_back(Vec3i::Constant(7));
    css.push_back(2.0f);
    css.push_back(2.0f);
    css.push_back(2.0f);
    css.push_back(2.0f);
    css.push_back(2.0f);
    css.push_back(2.0f);
    css.push_back(2.0f);
    css.push_back(2.0f);
    os.appendCells(tvI.begin(),tvI.end(),VTKWriter<scalarD>::POINT);

    tvI.clear();
    tvI.push_back(Vec3i(4,8,0));
    tvI.push_back(Vec3i(5,9,0));
    tvI.push_back(Vec3i(6,10,0));
    tvI.push_back(Vec3i(7,11,0));
    css.push_back(3.0f);
    css.push_back(3.0f);
    css.push_back(3.0f);
    css.push_back(3.0f);
    os.appendCells(tvI.begin(),tvI.end(),VTKWriter<scalarD>::LINE);
    os.appendCustomData("Color",css.begin(),css.end());
    _id++;
}

ClothCollision::NarrowNode::NarrowNode():Serializable(-1) {}
bool ClothCollision::NarrowNode::read(std::istream& is,IOData* dat)
{
    ASSERT_MSG(false,"Not Supported!");
    return false;
}
bool ClothCollision::NarrowNode::write(std::ostream& os,IOData* dat) const
{
    ASSERT_MSG(false,"Not Supported!");
    return false;
}
void ClothCollision::NarrowNode::buildBVH()
{
    sizeType nrT=(sizeType)_mesh->_tss.size();
    _bvh.resize(nrT);
    for(sizeType i=0; i<nrT; i++) {
        _bvh[i]._nrCell=1;
        _bvh[i]._parent=-1;
        _bvh[i]._cell=i;
        _mesh->_tss[i]->_activeTag=&_active;
    }
    refit(false);
    GEOM::BVHBuilder<GEOM::Node<sizeType,BBOX>,3>().buildBVH(_bvh);
    for(sizeType i=nrT; i<(sizeType)_bvh.size(); i++)
        _bvh[i]._cell=-1;
    //refit();
    //writeBVHByLevel(_bvh,boost::shared_ptr<ClothMesh::ClothTriangle>());
}
ClothCollision::BBOX ClothCollision::NarrowNode::refit(bool useLastPos)
{
    _vbb.resize(_mesh->_vss.size());
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<(sizeType)_mesh->_vss.size(); i++) {
        _vbb[i].reset();
        _vbb[i].setUnion(_mesh->_vss[i]->_pos.cast<scalar>());
        if(useLastPos)
            _vbb[i].setUnion(_mesh->_vss[i]->_lastPos.cast<scalar>());
        _vbb[i].enlarged((scalar)_thickness);
    }
    _ebb.resize(_mesh->_ess.size());
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<(sizeType)_mesh->_ess.size(); i++) {
        _ebb[i].reset();
        _ebb[i].setUnion(_vbb[_mesh->_ess[i]->_v[0]->_index]);
        _ebb[i].setUnion(_vbb[_mesh->_ess[i]->_v[1]->_index]);
    }
    for(sizeType i=0; i<(sizeType)_bvh.size(); i++) {
        _bvh[i]._bb.reset();
        if(_bvh[i]._cell >= 0) {
            const ClothMesh::ClothTriangle& t=*(_mesh->_tss[_bvh[i]._cell]);
            _bvh[i]._bb.setUnion(_ebb[t._e[0]->_index]);
            _bvh[i]._bb.setUnion(_ebb[t._e[1]->_index]);
            _bvh[i]._bb.setUnion(_ebb[t._e[2]->_index]);
        } else {
            _bvh[i]._bb.setUnion(_bvh[_bvh[i]._l]._bb);
            _bvh[i]._bb.setUnion(_bvh[_bvh[i]._r]._bb);
        }
    }
    if(!_active.empty()) {
        sizeType nrLeave=((sizeType)_bvh.size()+1)/2;
        ASSERT(_active.size() == nrLeave || _active.size() == _bvh.size());

        _active.resize(_bvh.size());
        for(sizeType i=nrLeave; i<(sizeType)_bvh.size(); i++)
            _active[i]=_active[_bvh[i]._l]||_active[_bvh[i]._r];
    }
    if(_bvh.empty())
        return BBOX();
    else return _bvh.back()._bb;
}
bool ClothCollision::NarrowNode::operator<(const NarrowNode& other) const
{
    return _mesh<other._mesh;
}

static FORCE_INLINE void adjustNormalSign(const Vec3d v[4],const Vec4d& omg,Vec3d& n)
{
    if(n.dot(v[0]*omg[0]+v[1]*omg[1]+
             v[2]*omg[2]+v[3]*omg[3]) < 0.0f)
        n*=-1.0f;
}
static FORCE_INLINE void calcCoef(const Vec3d pt[4],const Vec3d v[4],scalarD c[4])
{
    Vec3d x21=pt[1]-pt[0];
    Vec3d x31=pt[2]-pt[0];
    Vec3d x41=pt[3]-pt[0];

    Vec3d v21=v[1]-v[0];
    Vec3d v31=v[2]-v[0];
    Vec3d v41=v[3]-v[0];

    Vec3d x21Xx31=x21.cross(x31);
    Vec3d v21Xv31=v21.cross(v31);
    Vec3d v21Cx31_x21Cv31=v21.cross(x31)+x21.cross(v31);

    c[0]=x21Xx31.dot(x41);
    c[1]=x21Xx31.dot(v41)+v21Cx31_x21Cv31.dot(x41);
    c[2]=v21Xv31.dot(x41)+v21Cx31_x21Cv31.dot(v41);
    c[3]=v21Xv31.dot(v41);
}
static FORCE_INLINE scalarD eval(const scalarD c[4],scalarD s)
{
    scalarD ret=c[3]*s+c[2];
    ret=ret*s+c[1];
    ret=ret*s+c[0];
    return ret;
}
static FORCE_INLINE bool secant(scalarD l,scalarD r,const scalarD c[4],scalarD& s, scalarD tol)
{
    //out of range
    if(l>r)return false;

    //already close to zero
    scalarD el=eval(c,l);
    if(std::abs(el) < tol) {
        s=l;
        return true;
    }

    //already close to zero
    scalarD er=eval(c,r);
    if(std::abs(er) < tol) {
        s=r;
        return true;
    }

    //no zero in interval
    if(el*er > 0.0f)
        return false;

    //secant search
    scalarD m,em;
    for(int i=0; i<1000; i++) {
        m=(el*r-er*l)/(el-er);
        em=eval(c,m);
        if(std::abs(em) < tol || std::abs(r-l) < ClothCollision::_timeRes) {
            s=m;
            return true;
        }

        if(el*em < 0.0f) {
            r=m;
            er=em;
        } else {
            l=m;
            el=em;
        }
    }
    return false;
}
static FORCE_INLINE int solveCubicAllSecant(scalarD l,scalarD r,const scalarD c[4],scalarD s[4])
{
    scalarD tol=1E-7f*(std::abs(c[0])+std::abs(c[1])+std::abs(c[2])+std::abs(c[3]));
    //determine the minimum and maximum
    //c[1] + (2*c[2]) * x + (3*c[3]) * x^2
    //sol=(-c[2] (+/-) sqrt(c[2]*c[2]-3*c[3]*c[1]) )/(3*c[3])
    scalarD sol0,sol1;
    scalarD deter=sqrt(c[2]*c[2]-3.0f*c[3]*c[1]);
    if(deter < 0.0f) {
        return secant(l,r,c,s[0],tol) ? 1 : 0;
    } else {
        sol0=(-c[2]-deter)/(3.0f*c[3]);
        sol1=(-c[2]+deter)/(3.0f*c[3]);
        if(sol0>sol1)std::swap(sol0,sol1);
        //scant search in each segment
        int nrSol=0;
        if(secant(l,
                  std::min<scalarD>(sol0,r),
                  c,s[nrSol],tol))nrSol++;
        if(secant(std::max<scalarD>(l,sol0),
                  std::min<scalarD>(r,sol1),
                  c,s[nrSol],tol))nrSol++;
        if(secant(std::max<scalarD>(l,sol1),
                  r,
                  c,s[nrSol],tol))nrSol++;
        return nrSol;
    }
}
static FORCE_INLINE bool testVT(const Vec3d& v,const Vec3d& t0,const Vec3d& t1,const Vec3d& t2,Vec4d& omega,Vec3d& n)
{
    Vec3d t13=t0-t2;
    Vec3d t23=t1-t2;
    n=t13.cross(t23);

    scalarD area=n.norm();
    n/=area;

    if(abs((v-t0).dot(n)) > ClothCollision::_thickness)
        return false;
    else {
        scalarD delta=ClothCollision::_thickness/std::max(std::sqrt(area),ClothCollision::_rounding);
        Vec3d t43=v-t2;

        Vec2d RHS(t13.dot(t43),t23.dot(t43));
        Mat2d LHS;
        LHS(0,0)=t13.dot(t13);
        LHS(1,0)=LHS(0,1)=t13.dot(t23);
        LHS(1,1)=t23.dot(t23);

        omega.block<2,1>(0,0)=LHS.inverse()*RHS;
        omega[2]=1.0f-omega[0]-omega[1];
        return omega[0] >= -delta &&
               omega[1] >= -delta &&
               omega[2] >= -delta &&
               omega[0] <= 1.0f+delta &&
               omega[1] <= 1.0f+delta &&
               omega[2] <= 1.0f+delta;
    }
}
static FORCE_INLINE bool testEE(const Vec3d& eA0,const Vec3d& eA1,const Vec3d& eB0,const Vec3d& eB1,Vec4d& omega,Vec3d& n)
{
    Vec3d t21=eA1-eA0;
    Vec3d t43=eB1-eB0;

    scalarD delta=ClothCollision::_thickness/std::max(std::max(t21.norm(),t43.norm()),ClothCollision::_rounding);
    Vec3d nt21=t21/std::max(t21.norm(),ClothCollision::_rounding);
    Vec3d nt43=t43/std::max(t43.norm(),ClothCollision::_rounding);
    n=(nt21).cross(nt43);

    scalarD nNorm=n.norm();
    if(nNorm < ClothCollision::_rounding) {
        /*//ignore parallel edges
        //David Harmon pointed out that parallel edge constraint
        //can be detected by other point triangle test
        Vec3d dA=eA1-eA0;
        dA/=std::max(dA.norm(),_rounding);

        Vec3d n=eA0-eB0;
        n-=n.dot(dA)*dA;
        if(n.norm() > _thickness)
        	return false;

        scalarD IA[2]={eA0.dot(dA),eA1.dot(dA)};
        scalarD IB[2]={eB0.dot(dA),eB1.dot(dA)};

        scalarD SIB[2]={min(IB[0],IB[1]),max(IB[0],IB[1])};
        if(IA[0] >= SIB[0] && IA[0] <= SIB[1]){
        	omega[0]=1.0f;
        	omega[1]=0.0f;
        	omega[2]=(IA[0]-IB[1])/(IB[0]-IB[1]);
        	omega[3]=1.0f-omega[2];
        	return true;
        }else if(IA[1] >= SIB[0] && IA[1] <= SIB[1]){
        	omega[0]=0.0f;
        	omega[1]=1.0f;
        	omega[2]=(IA[1]-IB[1])/(IB[0]-IB[1]);
        	omega[3]=1.0f-omega[2];
        	return true;
        }

        scalarD SIA[2]={min(IA[0],IA[1]),max(IA[0],IA[1])};
        if(IB[0] >= SIA[0] && IB[0] <= SIA[1]){
        	omega[2]=1.0f;
        	omega[3]=0.0f;
        	omega[0]=(IB[0]-IA[1])/(IA[0]-IA[1]);
        	omega[1]=1.0f-omega[0];
        	return true;
        }else if(IB[1] >= SIA[0] && IB[1] <= SIA[1]){
        	omega[2]=0.0f;
        	omega[3]=1.0f;
        	omega[0]=(IB[1]-IA[1])/(IA[0]-IA[1]);
        	omega[1]=1.0f-omega[0];
        	return true;
        }*/
        return false;
    } else {
        n/=nNorm;
        Vec3d t31=eB0-eA0;

        Vec2d RHS(t21.dot(t31),-t43.dot(t31));
        Mat2d LHS;
        LHS(0,0)=t21.dot(t21);
        LHS(0,1)=-t21.dot(t43);
        LHS(1,0)=-t21.dot(t43);
        LHS(1,1)=t43.dot(t43);
        RHS=LHS.inverse()*RHS;
        if(RHS[0] >= -delta && RHS[1] >= -delta && RHS[0] <= 1.0f+delta && RHS[1] <= 1.0f+delta) {
            t21=eA0+t21*RHS[0];
            t43=eB0+t43*RHS[1];
            if((t21-t43).norm() < ClothCollision::_thickness) {
                omega[0]=1.0f-RHS[0];
                omega[1]=RHS[0];
                omega[2]=1.0f-RHS[1];
                omega[3]=RHS[1];
                return true;
            }
        }
        return false;
    }
}

bool testVT(boost::shared_ptr<ClothMesh::ClothVertex> V1,boost::shared_ptr<ClothMesh::ClothTriangle> T2,ClothCollision::CollisionHandler& handler)
{
    //geometry info
    const ClothMesh::ClothVertex& v0=*(T2->getV0());
    const ClothMesh::ClothVertex& v1=*(T2->getV1());
    const ClothMesh::ClothVertex& v2=*(T2->getV2());
    Vec3d n,pt[4]= {
        V1->_lastPos,
        v0._lastPos,
        v1._lastPos,
        v2._lastPos
    };
    Vec3d v[4]= {
        V1->_pos-V1->_lastPos,
        v0._pos-v0._lastPos,
        v1._pos-v1._lastPos,
        v2._pos-v2._lastPos,
    };

    //CCD test
    scalarD c[4];
    calcCoef(pt,v,c);

    scalarD s[4];
    scalarD cofl=0.0f,cofr=1.0f;
    int nr=solveCubicAllSecant(cofl,cofr,c,s);

    //check result
    Vec3d pos[4];
    Vec4d omega;
    for(int i = 0; i < nr; ++i) {
        for(int j = 0; j < 4; ++j)
            pos[j] = pt[j] + v[j]*s[i];
        if(testVT(pos[0], pos[1], pos[2], pos[3], omega, n)) {
            Vec4d omg(1.0f,-omega[0],-omega[1],-omega[2]);
            adjustNormalSign(v,-omg,n);
            handler.handle(V1,T2,n,omg,s[i]);
            return true;
        }
    }
    return false;
}
bool testEE(boost::shared_ptr<ClothMesh::ClothEdge> E1,boost::shared_ptr<ClothMesh::ClothEdge> E2,ClothCollision::CollisionHandler& handler)
{
    //geometry info
    const ClothMesh::ClothVertex& E10=*(E1->_v[0]);
    const ClothMesh::ClothVertex& E11=*(E1->_v[1]);
    const ClothMesh::ClothVertex& E20=*(E2->_v[0]);
    const ClothMesh::ClothVertex& E21=*(E2->_v[1]);
    Vec3d n,pt[4]= {
        E10._lastPos,
        E11._lastPos,
        E20._lastPos,
        E21._lastPos
    };
    Vec3d v[4]= {
        E10._pos-E10._lastPos,
        E11._pos-E11._lastPos,
        E20._pos-E20._lastPos,
        E21._pos-E21._lastPos
    };

    //CCD test
    scalarD c[4];
    calcCoef(pt,v,c);

    scalarD s[4];
    scalarD cofl = 0.0f, cofr = 1.0f;
    int nr = solveCubicAllSecant(cofl, cofr, c, s);

    //check result
    Vec3d pos[4];
    Vec4d omega;
    for(int i = 0; i < nr; ++i) {
        for(int j = 0; j < 4; ++j)
            pos[j] = pt[j] + v[j]*s[i];
        if(testEE(pos[0],pos[1],pos[2],pos[3],omega,n)) {
            Vec4d omg(omega[0],omega[1],-omega[2],-omega[3]);
            adjustNormalSign(v,-omg,n);
            handler.handle(E1,E2,n,omg,s[i]);
            return true;
        }
    }
    return false;
}
sizeType ClothCollision::nrMesh() const
{
    return ((sizeType)_bvh.size()+1)/2;
}
const ClothMesh& ClothCollision::getMesh(sizeType i) const
{
    return *(_bvh[i]._cell->_mesh);
}
void ClothCollision::updateMesh(boost::shared_ptr<ClothMesh> mesh)
{
    for(sizeType i=0; i<(sizeType)_bvh.size(); i++)
        if(_bvh[i]._cell && _bvh[i]._cell->_mesh == mesh)
            _bvh[i]._cell->buildBVH();
}
void ClothCollision::addMesh(boost::shared_ptr<ClothMesh> mesh)
{
    for(sizeType i=0; i<(sizeType)_bvh.size(); i++)
        if(_bvh[i]._cell && _bvh[i]._cell->_mesh == mesh)
            return;
    for(sizeType i=0; i<(sizeType)_bvh.size();)
        if(_bvh[i]._cell) {
            _bvh[i]._nrCell=1;
            _bvh[i]._parent=-1;
            _bvh[i]._bb=_bvh[i]._cell->refit();
            i++;
        } else {
            _bvh[i]=_bvh.back();
            _bvh.pop_back();
        }
    _bvh.push_back(GEOM::Node<boost::shared_ptr<NarrowNode>,BBOX>());
    _bvh.back()._nrCell=1;
    _bvh.back()._parent=-1;
    _bvh.back()._cell.reset(new NarrowNode);
    _bvh.back()._cell->_mesh=mesh;
    _bvh.back()._cell->buildBVH();
    _bvh.back()._bb=_bvh.back()._cell->refit();
    //parity check
    for(sizeType i=0; i<(sizeType)_bvh.size(); i++) {
        if(_bvh[i]._cell->_mesh->_vss.empty())
            continue;
        if(i == 0) {
            ASSERT_MSG(_bvh[i]._cell->_mesh->_vss[0]->_type == ClothMesh::CLOTH_MESH,"We only allow one ClothMesh be inserted first!");
        } else {
            ASSERT_MSG(_bvh[i]._cell->_mesh->_vss[0]->_type != ClothMesh::CLOTH_MESH,"We only allow one ClothMesh be inserted first!");
        }
    }
    //rebuild BVH
    GEOM::BVHBuilder<GEOM::Node<boost::shared_ptr<NarrowNode>,BBOX>,3>().buildBVH(_bvh);
}
void ClothCollision::delMesh(boost::shared_ptr<ClothMesh> mesh)
{
    while(!_bvh.empty() && !_bvh.back()._cell)
        _bvh.pop_back();
    for(sizeType i=0; i<(sizeType)_bvh.size();)
        if(_bvh[i]._cell->_mesh == mesh) {
            _bvh[i]=_bvh.back();
            _bvh.pop_back();
        } else {
            _bvh[i]._nrCell=1;
            _bvh[i]._parent=-1;
            _bvh[i]._cell->refit();
            i++;
        }
    GEOM::BVHBuilder<GEOM::Node<boost::shared_ptr<NarrowNode>,BBOX>,3>().buildBVH(_bvh);
}
void ClothCollision::collide(CollisionHandler& handler,bool useActive)
{
    sizeType nrFound,nrFiltered,nrReal;
    //refit all
    for(sizeType i=0; i<(sizeType)_bvh.size(); i++)
        if(_bvh[i]._cell)
            _bvh[i]._bb=_bvh[i]._cell->refit();
        else {
            _bvh[i]._bb.reset();
            _bvh[i]._bb.setUnion(_bvh[_bvh[i]._l]._bb);
            _bvh[i]._bb.setUnion(_bvh[_bvh[i]._r]._bb);
        }

    //BroadPhase
    _cache.clear();
    GEOM::BVHQuery<boost::shared_ptr<NarrowNode>,BBOX> q(_bvh,3,boost::shared_ptr<NarrowNode>());
    q.interBodyQuery(q,*this);
    std::sort(_cache.begin(),_cache.end());
    _cache.erase(std::unique(_cache.begin(),_cache.end()),_cache.end());

    //NarrowPhase
    _cacheVT.clear();
    _cacheEE.clear();
    for(sizeType i=0; i<(sizeType)_cache.size(); i++) {
        _activeCache=_cache[i];
//#define BRUTE_FORCE
#ifdef BRUTE_FORCE
        vector<boost::shared_ptr<ClothMesh::ClothTriangle> >& tssA=_activeCache._A->_mesh->_tss;
        vector<boost::shared_ptr<ClothMesh::ClothTriangle> >& tssB=_activeCache._B->_mesh->_tss;
        for(sizeType r=0; r<(sizeType)tssA.size(); r++) {
            BBOX bbR;
            bbR.setUnion(tssA[r]->getV0()->_pos.cast<scalar>());
            bbR.setUnion(tssA[r]->getV1()->_pos.cast<scalar>());
            bbR.setUnion(tssA[r]->getV2()->_pos.cast<scalar>());
            bbR.setUnion(tssA[r]->getV0()->_lastPos.cast<scalar>());
            bbR.setUnion(tssA[r]->getV1()->_lastPos.cast<scalar>());
            bbR.setUnion(tssA[r]->getV2()->_lastPos.cast<scalar>());
            for(sizeType c=0; c<(sizeType)tssB.size(); c++) {
                BBOX bbC;
                bbC.setUnion(tssB[c]->getV0()->_pos.cast<scalar>());
                bbC.setUnion(tssB[c]->getV1()->_pos.cast<scalar>());
                bbC.setUnion(tssB[c]->getV2()->_pos.cast<scalar>());
                bbC.setUnion(tssB[c]->getV0()->_lastPos.cast<scalar>());
                bbC.setUnion(tssB[c]->getV1()->_lastPos.cast<scalar>());
                bbC.setUnion(tssB[c]->getV2()->_lastPos.cast<scalar>());
                if(bbR.intersect(bbC))
                    onCell(_cache[i]._A->_bvh[r],_cache[i]._B->_bvh[c]);
            }
        }
#else
        GEOM::BVHQuery<sizeType,BBOX> q1(_cache[i]._A->_bvh,3,-1);
        GEOM::BVHQuery<sizeType,BBOX> q2(_cache[i]._B->_bvh,3,-1);
        if(useActive && !_cache[i]._A->_active.empty())
            q1._active=&(_cache[i]._A->_active);
        if(useActive && !_cache[i]._B->_active.empty())
            q2._active=&(_cache[i]._B->_active);
        q1.interBodyQuery(q2,*this);
#endif
    }

    //Fine-Grained TV Test
    nrFound=_cacheVT.size();
    std::sort(_cacheVT.begin(),_cacheVT.end());
    _cacheVT.erase(std::unique(_cacheVT.begin(),_cacheVT.end()),_cacheVT.end());
    nrFiltered=_cacheVT.size();
    nrReal=0;
    for(sizeType i=0; i<(sizeType)_cacheVT.size(); i++)
        if(testVT(_cacheVT[i]._A,_cacheVT[i]._B,handler))
            nrReal++;
    INFOV("VT Collision: (%d,%d,%d)",nrFound,nrFiltered,nrReal);

    //Fine-Grained EE Test
    nrFound=_cacheEE.size();
    std::sort(_cacheEE.begin(),_cacheEE.end());
    _cacheEE.erase(std::unique(_cacheEE.begin(),_cacheEE.end()),_cacheEE.end());
    nrFiltered=_cacheEE.size();
    nrReal=0;
    for(sizeType i=0; i<(sizeType)_cacheEE.size(); i++)
        if(testEE(_cacheEE[i]._A,_cacheEE[i]._B,handler))
            nrReal++;
    INFOV("EE Collision: (%d,%d,%d)",nrFound,nrFiltered,nrReal);
}
void ClothCollision::restartActive()
{
    for(sizeType i=0; i<(sizeType)_bvh.size(); i++)
        if(_bvh[i]._cell) {
            NarrowNode& n=*(_bvh[i]._cell);
            n._active.assign(((sizeType)n._bvh.size()+1)/2,false);
        }
}
void ClothCollision::activate(const ClothMesh::ClothVertex* t)
{
    for(sizeType i=0; i<(sizeType)t->_oneRing.size(); i++) {
        ClothMesh::ClothTriangle& neigh=*(t->_oneRing[i]);
        vector<bool>& tag=*(neigh._activeTag);
        if(!tag.empty())
            tag[neigh._index]=true;
    }
}
void ClothCollision::onCell(const GEOM::Node<boost::shared_ptr<NarrowNode>,BBOX>& A,const GEOM::Node<boost::shared_ptr<NarrowNode>,BBOX>& B)
{
    bool rigidA=A._cell->_mesh->_vss[0]->_type&ClothMesh::RIGID_MESH;
    bool rigidB=B._cell->_mesh->_vss[0]->_type&ClothMesh::RIGID_MESH;
    if(!rigidA || !rigidB)
        _cache.push_back(Cache<NarrowNode,NarrowNode>(A._cell,B._cell));
}
static FORCE_INLINE void addColl
(boost::shared_ptr<ClothMesh::ClothEdge> eA,boost::shared_ptr<ClothMesh::ClothEdge> eB,
 vector<ClothCollision::Cache<ClothMesh::ClothEdge,ClothMesh::ClothEdge> >& cacheEE,bool sameMesh)
{
    if(!sameMesh ||
            (eA->_v[0]->_index != eB->_v[0]->_index &&
             eA->_v[0]->_index != eB->_v[1]->_index &&
             eA->_v[1]->_index != eB->_v[0]->_index &&
             eA->_v[1]->_index != eB->_v[1]->_index))
        cacheEE.push_back(ClothCollision::Cache<ClothMesh::ClothEdge,ClothMesh::ClothEdge>(eA,eB));
}
static FORCE_INLINE void addColl
(boost::shared_ptr<ClothMesh::ClothVertex> eA,boost::shared_ptr<ClothMesh::ClothTriangle> eB,
 vector<ClothCollision::Cache<ClothMesh::ClothVertex,ClothMesh::ClothTriangle> >& cacheVT,bool sameMesh)
{
    if(!sameMesh ||
            (eA->_index != eB->getV0()->_index &&
             eA->_index != eB->getV1()->_index &&
             eA->_index != eB->getV2()->_index))
        cacheVT.push_back(ClothCollision::Cache<ClothMesh::ClothVertex,ClothMesh::ClothTriangle>(eA,eB));
}
void ClothCollision::onCell(const GEOM::Node<sizeType,BBOX>& nA,const GEOM::Node<sizeType,BBOX>& nB)
{
    const NarrowNode& nodeA=*(_activeCache._A);
    const NarrowNode& nodeB=*(_activeCache._B);
    boost::shared_ptr<ClothMesh::ClothTriangle> tA=nodeA._mesh->_tss[nA._cell];
    boost::shared_ptr<ClothMesh::ClothTriangle> tB=nodeB._mesh->_tss[nB._cell];

    //push EE
    boost::shared_ptr<ClothMesh::ClothEdge> AE[3]= {tA->_e[0],tA->_e[1],tA->_e[2]};
    boost::shared_ptr<ClothMesh::ClothEdge> BE[3]= {tB->_e[0],tB->_e[1],tB->_e[2]};
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++) {
            //ASSERT(nodeA._mesh->_ess[AE[i]->_index] == AE[i]);
            //ASSERT(nodeB._mesh->_ess[BE[i]->_index] == BE[i]);
            if(AE[i] != BE[j] && nodeA._ebb[AE[i]->_index].intersect(nodeB._ebb[BE[j]->_index]))
                addColl(AE[i],BE[j],_cacheEE,nodeA._mesh == nodeB._mesh);
        }

    //push TB
    boost::shared_ptr<ClothMesh::ClothVertex> AV[3]= {tA->getV0(),tA->getV1(),tA->getV2()};
    boost::shared_ptr<ClothMesh::ClothVertex> BV[3]= {tB->getV0(),tB->getV1(),tB->getV2()};
    for(int i=0; i<3; i++) {
        //ASSERT(nodeA._mesh->_vss[AV[i]->_index] == AV[i]);
        //ASSERT(nodeB._mesh->_vss[BV[i]->_index] == BV[i]);
        if(nodeA._vbb[AV[i]->_index].intersect(nB._bb))addColl(AV[i],tB,_cacheVT,nodeA._mesh == nodeB._mesh);
        if(nodeB._vbb[BV[i]->_index].intersect(nA._bb))addColl(BV[i],tA,_cacheVT,nodeA._mesh == nodeB._mesh);
    }
}
scalarD ClothCollision::_thickness=1E-3f;
scalarD ClothCollision::_rounding=1E-6f;
scalarD ClothCollision::_timeRes=1E-4f;
