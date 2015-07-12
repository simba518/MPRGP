#include "FEMCollision.h"
#include "FEMCollider.h"
#include "FEMMesh.h"
#include "FEMGeom.h"
#include "FEMBVHBuilder.h"
#include "FEMUtils.h"
#include "Heap.h"
#include <stack>
#include "CollisionDetection.h"

USE_PRJ_NAMESPACE

void collideCC
(vector<FEMCache<FEMCell,FEMVertex> >& cache,
 boost::shared_ptr<FEMBody> bA,boost::shared_ptr<FEMCell> nA,
 boost::shared_ptr<FEMBody> bB,boost::shared_ptr<FEMCell> nB,
 const BBox<scalar>& bbA,bool surface,sizeType dim)
{
    //early out
    if(nA == nB)
        return;

    //we don't test connected tets
    for(sizeType i=0; i<dim+1; i++)
        for(sizeType j=0; j<dim+1; j++)
            if(nA->_v[i] == nB->_v[j])
                return;

    for(sizeType j=0; j<dim+1; j++) {
        //not contain
        if(!bbA.contain(nB->_v[j]->_pos,dim))
            continue;
        //not surface node
        if(surface && !nB->_v[j]->_surface)
            continue;
        cache.push_back(FEMCache<FEMCell,FEMVertex>(bA,nA,bB,nB->_v[j]));
    }
}
//Default Collision Solver //slow but simple implementation
class DefaultFEMBody : public FEMBody
{
public:
    DefaultFEMBody():FEMBody() {}
    virtual FEMBody& operator=(const FEMBody& other) {
        const DefaultFEMBody* otherDef=static_cast<const DefaultFEMBody*>(&other);
        ASSERT_MSG(otherDef,"Different Body Type!");
        FEMBody::operator=(other);
        return *this;
    }
    boost::shared_ptr<Serializable> copy() const {
        return boost::shared_ptr<Serializable>(new DefaultFEMBody);
    }
};
sizeType FEMCollision::dim() const
{
    return _mesh->dim();
}
boost::shared_ptr<FEMBody> FEMCollision::createBody() const
{
    return boost::shared_ptr<FEMBody>(new DefaultFEMBody);
}
boost::shared_ptr<FEMCollision> FEMCollision::copy() const
{
    return boost::shared_ptr<FEMCollision>(new FEMCollision);
}
void FEMCollision::updateMesh() {}
void FEMCollision::rayCastPSet(const LineSegTpl<scalar>& l,scalar rad,FEMInterp& cd,scalar& dist) const
{
    sizeType nrB=_mesh->nrB();
    for(sizeType b=0; b<nrB; b++) {
        const FEMBody& body=_mesh->getB(b);
        const ParticleSetN& pSet=body.getPSet();
        sizeType nrP=pSet.size();
        for(sizeType i=0; i<nrP; i++) {
            Vec3 cp,b,dir=(l._y-l._x).normalized();
            scalar sqrDist;
            l.calcPointDist(pSet[i]._pos,sqrDist,cp,b);
            if(std::sqrt(sqrDist) < rad) {
                scalar d=(pSet[i]._pos-l._x).dot(dir);
                if(d < dist) {
                    dist=d;
                    cd=body.getI(i);
                }
            }
        }
    }
}
void FEMCollision::rayCastMesh(const LineSegTpl<scalar>& l,FEMInterp& cd,scalar& dist) const
{
    sizeType nrB=_mesh->nrB();
    for(sizeType b=0; b<nrB; b++) {
        const FEMBody& body=_mesh->getB(b);
        sizeType nrC=body.nrC();
        for(sizeType i=0; i<nrC; i++) {
            const FEMCell& c=body.getC(i);
            Vec4 bary;
            Vec2 lb;
            bool hasColl;
            if(_mesh->dim() == 3) {
                TetrahedronTpl<scalar> tet(c[0],c[1],c[2],c[3]);
                hasColl=tet.calcLineDist(l,bary,lb);
            } else {
                Vec3 bary3;
                TriangleTpl<scalar> tri(c[0],c[1],c[2]);
                hasColl=tri.calcLineDist(l,bary3,lb);
                bary.block<3,1>(0,0)=bary3;
            }
            if(hasColl) {
                scalar d=(l._x*lb[0]+l._y*lb[1]-l._x).norm();
                if(d < dist) {
                    dist=d;
                    cd._coef=bary;
                    cd._cell=body.getCPtr(i);
                }
            }
        }
    }
}
void FEMCollision::collideMesh(FEMCollider& coll,bool surface)
{
    _cache.clear();
    sizeType nrB=_mesh->nrB();
    for(sizeType bi=0; bi<nrB; bi++)
        for(sizeType bj=0; bj<nrB; bj++) {
            boost::shared_ptr<FEMBody> bodyA=_mesh->getBPtr(bi);
            boost::shared_ptr<FEMBody> bodyB=_mesh->getBPtr(bj);
            for(sizeType i=0; i<bodyA->nrC(); i++)
                for(sizeType j=0; j<bodyB->nrC(); j++) {
                    BBox<scalar> bbA=bodyA->getC(i).getBB();
                    BBox<scalar> bbB=bodyB->getC(j).getBB();
                    if(!bbA.intersect(bbB,_mesh->dim()))continue;
                    collideCC(_cache,
                              bodyA,bodyA->getCPtr(i),
                              bodyB,bodyB->getCPtr(j),
                              bbA,surface,_mesh->dim());
                }
        }

    std::sort(_cache.begin(),_cache.end());
    sizeType nrUnique=0;
    sizeType nrColl=(sizeType)_cache.size();
    for(sizeType i=0; i<nrColl; i++)
        if(i == 0 || _cache[i] != _cache[i-1]) {
            const FEMCache<FEMCell,FEMVertex>& c=_cache[i];
            FEMInterp I;
            if(c._A->contain(c._B->_pos,I))
                coll.handle(c._bA,c._A,c._bB,c._B,I._coef);
            nrUnique++;
        }
    //INFOV("Nr Cache: %lu, NrUnique: %lu",_cache.size(),nrUnique)
}
void FEMCollision::collideGeom(const FEMGeom& other,FEMCollider& coll,bool surface)
{
    sizeType nrB=_mesh->nrB();
    for(sizeType b=0; b<nrB; b++) {
        boost::shared_ptr<FEMBody> body=_mesh->getBPtr(b);
        sizeType nrV=body->nrV();
        sizeType nrG=other.nrG();
        for(sizeType i=0; i<nrV; i++) {
            if(surface && !body->getV(i)._surface)continue;
            Vec3 n,pt=body->getV(i)._pos;
            for(sizeType j=0; j<nrG; j++)
                if(other.getG(j).dist(pt,n))
                    coll.handle(body,body->getVPtr(i),n);
        }
    }
}

//BVHCollisionBody
class VBVHFEMBody : public FEMBody
{
protected:
    struct Face {
        Face() {}
        Face(const Vec3i& vert,const Vec2i& cell):_vert(vert),_cell(cell) {}
        scalar cost(const vector<Node<sizeType> >& bvh,sizeType& t0,sizeType& t1,sizeType dim) const {
            t0=_cell[0];
            while(bvh[t0]._parent != -1)t0=bvh[t0]._parent;
            t1=_cell[1];
            while(bvh[t1]._parent != -1)t1=bvh[t1]._parent;
            if(t0 == t1)
                return -1.0f;

            /*scalar retL,retR,ret;
            if(dim == 2){
            	retL=SurfaceArea<2>::area(bvh[t0]._bb);
            	retR=SurfaceArea<2>::area(bvh[t1]._bb);
            	ret=SurfaceArea<2>::area(bvh[t0]._bb.getUnion(bvh[t1]._bb));
            }else{
            	retL=SurfaceArea<3>::area(bvh[t0]._bb);
            	retR=SurfaceArea<3>::area(bvh[t1]._bb);
            	ret=SurfaceArea<3>::area(bvh[t0]._bb.getUnion(bvh[t1]._bb));
            }ret/=(retL+retR);
            return ret;*/

            scalar ret=(scalar)(bvh[t0]._nrCell+bvh[t1]._nrCell);
            Vec3 ext=bvh[t0]._bb.getUnion(bvh[t1]._bb).getExtent();
            for(int i=0; i<dim; i++)ret*=ext[i];
            return ret;
        }
        Vec3i _vert;
        Vec2i _cell;
    };
    void mergeNode(const vector<Face>& faces) {
        scalar res;
        sizeType t0,t1;
        vector<scalarD> result;
        vector<sizeType> heap,heapOff;
        sizeType nrF=(sizeType)faces.size();
        result.assign(nrF,0.0f);
        heapOff.assign(nrF,-1);
        for(sizeType i=0; i<nrF; i++) {
            result[i]=faces[i].cost(_bvh,t0,t1,dim());
            pushHeapDef(result,heapOff,heap,i);
        }

        //merge using face
        while(heap.size() > 1) {
            while(heap.size() > 1) {
                res=faces[heap[0]].cost(_bvh,t0,t1,dim());
                if(res < 0.0f) {
                    popHeapDef(result,heapOff,heap,-1);
                } else if(res != result[heap[0]]) {
                    t0=popHeapDef(result,heapOff,heap,-1);
                    result[t0]=res;
                    pushHeapDef(result,heapOff,heap,t0);
                } else {
                    popHeapDef(result,heapOff,heap,-1);
                    break;
                }
            }
            if(res >= 0.0f && t0 != t1) {
                _bvh[t0]._parent=(sizeType)_bvh.size();
                _bvh[t1]._parent=(sizeType)_bvh.size();

                Node<sizeType> n;
                n._cell=verbose();
                n._bb=_bvh[t0]._bb;
                n._bb.setUnion(_bvh[t1]._bb);
                n._l=t0;
                n._r=t1;
                n._parent=-1;
                n._nrCell=_bvh[n._l]._nrCell+_bvh[n._r]._nrCell;
                _bvh.push_back(n);
            }
        }
    }
    struct DistCallback {
        DistCallback(const FEMBody& body,const Vec3& pt,sizeType dim):_body(body),_pt(pt),_dim(dim) {}
        void updateDist(const Node<sizeType>& node,const Vec3& pt,Vec3& cp,Vec3& n,scalar& dist) const {
            if(_dim == 2)updateDist2D(_body.getC(node._cell),_pt,dist);
            else updateDist3D(_body.getC(node._cell),_pt,dist);
        }
        void updateDist2D(const FEMCell& c,const Vec3& pt,scalar& distNew) const {
            static const sizeType lid[3][2]= {{1,2},{0,2},{0,1}};
            Vec3 cp,b;
            scalar sqrDist;
            for(sizeType l=0; l<3; l++) {
                bool interior=false;
                for(sizeType lv=0; lv<2; lv++)
                    if(!c._v[lid[l][lv]]->_surface) {
                        interior=true;
                        break;
                    }
                if(interior)continue;

                LineSeg seg(c[lid[l][0]],c[lid[l][1]]);
                seg.calcPointDist(pt,sqrDist,cp,b);
                distNew=std::min<scalar>(sqrt(sqrDist),distNew);
            }
        }
        void updateDist3D(const FEMCell& c,const Vec3& pt,scalar& distNew) const {
            static const sizeType lid[4][3]= {{1,2,3},{0,2,3},{0,1,3},{0,1,2}};
            Vec3 cp,b;
            scalar sqrDist;
            for(sizeType l=0; l<4; l++) {
                bool interior=false;
                for(sizeType lv=0; lv<3; lv++)
                    if(!c._v[lid[l][lv]]->_surface) {
                        interior=true;
                        break;
                    }
                if(interior)continue;

                Triangle seg(c[lid[l][0]],c[lid[l][1]],c[lid[l][2]]);
                seg.calcPointDist(pt,sqrDist,cp,b);
                distNew=std::min<scalar>(sqrt(sqrDist),distNew);
            }
        }
        scalar depth() const {
            return numeric_limits<scalar>::max();
        }
        const FEMBody& _body;
        const Vec3& _pt;
        sizeType _dim;
    };
public:
    VBVHFEMBody():FEMBody() {}
    virtual FEMBody& operator=(const FEMBody& other) {
        const VBVHFEMBody* otherBVH=static_cast<const VBVHFEMBody*>(&other);
        ASSERT_MSG(otherBVH,"Different Body Type!");
        FEMBody::operator=(other);
        _bvh=otherBVH->_bvh;
        return *this;
    }
    bool read(std::istream& is,IOData* dat) {
        FEMBody::read(is,dat);
        readVector(_bvh,is);
        return is.good();
    }
    bool write(std::ostream& os,IOData* dat) const {
        FEMBody::write(os,dat);
        writeVector(_bvh,os);
        return os.good();
    }
    boost::shared_ptr<Serializable> copy() const {
        return boost::shared_ptr<Serializable>(new VBVHFEMBody);
    }
    virtual void assemble() {
        FEMBody::assemble();
        buildBVH();
    }
    //collision
    virtual void rayCastPSet(const LineSegTpl<scalar>& l,scalar rad,FEMInterp& cd,scalar& dist) const {
        std::stack<sizeType> ss;
        ss.push(_bvh.size()-1);
        while(!ss.empty()) {
            const Node<sizeType>& node=_bvh[ss.top()];
            ss.pop();
            const BBox<scalar> bb=node._bb.enlarge(rad,dim());
            //exit if box is away from line
            scalar s,t;
            if(!bb.intersect(l._x,l._y,s,t,dim()))continue;
            if(l.length()*s > dist)continue;

            //test
            if(node._cell != verbose()) {
                boost::shared_ptr<FEMInterp> it=getCell(node)->_pSet;
                while(it) {
                    Vec3 cp,b,dir=(l._y-l._x).normalized();
                    scalar sqrDist;
                    l.calcPointDist(_pSet[it->_id]._pos,sqrDist,cp,b);
                    if(std::sqrt(sqrDist) < rad) {
                        scalar d=(_pSet[it->_id]._pos-l._x).dot(dir);
                        if(d < dist) {
                            dist=d;
                            cd=*it;
                        }
                    }
                    it=it->_next;
                }
            } else {
                ss.push(node._l);
                ss.push(node._r);
            }
        }
    }
    virtual void rayCastMesh(const LineSegTpl<scalar>& l,FEMInterp& cd,scalar& dist) const {
        std::stack<sizeType> ss;
        ss.push(_bvh.size()-1);
        while(!ss.empty()) {
            const Node<sizeType>& node=_bvh[ss.top()];
            ss.pop();
            //const scalar len=(bb._minC-bb._maxC).norm()*0.5f;
            //exit if box is away from line
            scalar s,t;
            if(!node._bb.intersect(l._x,l._y,s,t,dim()))continue;
            if(l.length()*s > dist)continue;

            //test
            if(node._cell != verbose()) {
                boost::shared_ptr<FEMCell> c=getCell(node);
                Vec4 bary;
                Vec2 lb;
                bool hasColl;
                if(dim() == 3) {
                    Tetrahedron tet((*c)[0],(*c)[1],(*c)[2],(*c)[3]);
                    hasColl=tet.calcLineDist(l,bary,lb);
                } else {
                    Vec3 bary3;
                    Triangle tri((*c)[0],(*c)[1],(*c)[2]);
                    hasColl=tri.calcLineDist(l,bary3,lb);
                    bary.block<3,1>(0,0)=bary3;
                }
                if(hasColl) {
                    scalar d=(l._x*lb[0]+l._y*lb[1]-l._x).norm();
                    if(d < dist) {
                        dist=d;
                        cd._coef=bary;
                        cd._cell=c;
                    }
                }
            } else {
                ss.push(node._l);
                ss.push(node._r);
            }
        }
    }
    //build BVH
    virtual void buildBVH() {
        //build leaves
        _bvh.resize(_css.size());
        for(sizeType i=0; i<(sizeType)_css.size(); i++) {
            _bvh[i]=Node<sizeType>();
            _bvh[i]._nrCell=1;
            _bvh[i]._cell=i;
            _bvh[i]._bb=_css[i]->getBB();
        }

        //build faces
        vector<std::pair<Vec3i,Vec2i> > faces;
        getFace(NULL,&faces);
        vector<Face> nodes((sizeType)faces.size());
        for(sizeType i=0; i<(sizeType)faces.size(); i++)
            nodes[i]=Face(faces[i].first,faces[i].second);

        //build BVH from bottom up
        sizeType sz=_bvh.size();
        mergeNode(nodes);
        // writeBVHByLevel<sizeType>(_bvh,verbose());
        ASSERT(_bvh.size() == sz*2-1);

        //build distance field
        sizeType nrV=(sizeType)_vss.size();
        BVHQuery<sizeType> query(_bvh,dim(),verbose());
        OMP_PARALLEL_FOR_
        for(sizeType i=0; i<nrV; i++) {
            FEMVertex& v=*(_vss[i]);
            if(v._surface)v._matDist=0.0f;
            else {
                Vec3 cp,n;
                v._matDist=numeric_limits<scalar>::max();
                DistCallback cb(*this,v._pos,dim());
                query.pointDistQuery(v._pos,cb,cp,n,v._matDist);
                v._matDist=-v._matDist;
            }
        }
    }
    virtual BBox<scalar> updateBVH() {
        sizeType nrC=(sizeType)_css.size();
        OMP_PARALLEL_FOR_
        for(sizeType i=0; i<nrC; i++)
            _bvh[i]._bb=_css[i]->getBB();
        BVHQuery<sizeType>(_bvh,dim(),verbose()).updateBVH();
        return _bvh.back()._bb;
    }
    virtual boost::shared_ptr<FEMCell> getCell(const Node<sizeType>& node) const {
        return _css[node._cell];
    }
    static FORCE_INLINE sizeType verbose() {
        return -1;
    }
    vector<Node<sizeType> > _bvh;
};
struct GeomCallback {
    GeomCallback(boost::shared_ptr<VBVHFEMBody> bA,boost::shared_ptr<FEMGeomCell> B,vector<FEMCache<FEMVertex,FEMGeomCell> >& cacheG,sizeType dim,bool surface)
        :_bA(bA),_B(B),_bbB(B->getBB()),_cacheG(cacheG),_dim(dim),_surface(surface) {}
    bool validNode(const Node<sizeType>& nA) const {
        return nA._bb.intersect(_bbB,_dim);
    }
    void updateDist(const Node<sizeType>& nA) const {
        boost::shared_ptr<FEMCell> cA=_bA->getCell(nA);
        for(sizeType v=0; v<4; v++) {
            if(!cA->_v[v])continue;
            if(_surface && !cA->_v[v]->_surface)continue;
            if(_bbB.contain(cA->_v[v]->_pos,_dim))
                _cacheG.push_back(FEMCache<FEMVertex,FEMGeomCell>(_bA,cA->_v[v],_bA,_B));
        }
    }
    boost::shared_ptr<VBVHFEMBody> _bA;
    boost::shared_ptr<FEMGeomCell> _B;
    vector<FEMCache<FEMVertex,FEMGeomCell> >& _cacheG;
    const BBox<scalar> _bbB;
    sizeType _dim;
    bool _surface;
};
struct MeshCallback {
    MeshCallback(boost::shared_ptr<FEMBody> bA,boost::shared_ptr<FEMBody> bB,vector<FEMCache<FEMCell,FEMVertex> >& cache,sizeType dim,bool surface)
        :_bA(bA),_bB(bB),_cache(cache),_dim(dim),_surface(surface) {}
    void onCell(const Node<sizeType>& nA,const Node<sizeType>& nB) {
        collideCC(_cache,_bA,_bA->getCPtr(nA._cell),_bB,_bB->getCPtr(nB._cell),nA._bb,_surface,_dim);
    }
    boost::shared_ptr<FEMBody> _bA,_bB;
    vector<FEMCache<FEMCell,FEMVertex> >& _cache;
    sizeType _dim;
    bool _surface;
};
boost::shared_ptr<FEMBody> BVHFEMCollision::createBody() const
{
    return boost::shared_ptr<FEMBody>(new VBVHFEMBody);
}
boost::shared_ptr<FEMCollision> BVHFEMCollision::copy() const
{
    return boost::shared_ptr<FEMCollision>(new BVHFEMCollision);
}
void BVHFEMCollision::updateMesh()
{
    sizeType nrB=_mesh->nrB();
    _bvh.reset(new BVH(nrB));
    for(sizeType i=0; i<nrB; i++) {
        Node<boost::shared_ptr<FEMBody> > n;
        n._l=n._r=n._parent=-1;
        n._nrCell=1;
        n._cell=_mesh->getBPtr(i);
        n._bb=boost::dynamic_pointer_cast<VBVHFEMBody>(n._cell)->updateBVH();
        (*_bvh)[i]=n;
    }
    buildBVH(*_bvh,_mesh->dim(),boost::shared_ptr<FEMBody>());
}
void BVHFEMCollision::rayCastPSet(const LineSegTpl<scalar>& l,scalar rad,FEMInterp& cd,scalar& dist) const
{
    const vector<Node<boost::shared_ptr<FEMBody> > >& bvh=*_bvh;
    std::stack<sizeType> ss;
    ss.push(bvh.size()-1);
    while(!ss.empty()) {
        const Node<boost::shared_ptr<FEMBody> >& n=bvh[ss.top()];
        ss.pop();
        const BBox<scalar> bb=n._bb;
        //const scalar len=(bb._minC-bb._maxC).norm()*0.5f;
        //exit if box is away from line
        scalar s,t;
        if(!bb.intersect(l._x,l._y,s,t,dim()))continue;
        if(l.length()*s > dist)continue;

        //test
        if(n._cell) {
            boost::dynamic_pointer_cast<VBVHFEMBody>(n._cell)->rayCastPSet(l,rad,cd,dist);
        } else {
            ss.push(n._l);
            ss.push(n._r);
        }
    }
}
void BVHFEMCollision::rayCastMesh(const LineSegTpl<scalar>& l,FEMInterp& cd,scalar& dist) const
{
    const vector<Node<boost::shared_ptr<FEMBody> > >& bvh=*_bvh;
    std::stack<sizeType> ss;
    ss.push(bvh.size()-1);
    while(!ss.empty()) {
        const Node<boost::shared_ptr<FEMBody> >& n=bvh[ss.top()];
        ss.pop();
        const BBox<scalar> bb=n._bb;
        //const scalar len=(bb._minC-bb._maxC).norm()*0.5f;
        //exit if box is away from line
        scalar s,t;
        if(!bb.intersect(l._x,l._y,s,t,dim()))continue;
        if(l.length()*s > dist)continue;

        //test
        if(n._cell) {
            boost::dynamic_pointer_cast<VBVHFEMBody>(n._cell)->rayCastMesh(l,cd,dist);
        } else {
            ss.push(n._l);
            ss.push(n._r);
        }
    }
}
void BVHFEMCollision::collideMesh(FEMCollider& coll,bool surface)
{
    vector<FEMCache<FEMBody,FEMBody> > cache;
    BVHQuery<boost::shared_ptr<FEMBody> > query(*_bvh,dim(),boost::shared_ptr<FEMBody>());
    query.broadphaseQuery(query,cache);

    _cache.clear();
    std::sort(cache.begin(),cache.end());
    for(sizeType i=0; i<(sizeType)cache.size(); i++)
        if(i == 0 || cache[i] != cache[i-1]) {
            boost::shared_ptr<VBVHFEMBody> bodyA=boost::dynamic_pointer_cast<VBVHFEMBody>(cache[i]._A);
            boost::shared_ptr<VBVHFEMBody> bodyB=boost::dynamic_pointer_cast<VBVHFEMBody>(cache[i]._B);
            BVHQuery<sizeType> query(bodyA->_bvh,dim(),bodyA->verbose());
            BVHQuery<sizeType> query2(bodyB->_bvh,dim(),bodyB->verbose());
            MeshCallback cb(cache[i]._A,cache[i]._B,_cache,dim(),surface);
            query.interBodyQuery(query2,cb);
        }

    std::sort(_cache.begin(),_cache.end());
    sizeType nrUnique=0;
    sizeType nrColl=(sizeType)_cache.size();
    for(sizeType i=0; i<nrColl; i++)
        if(i == 0 || _cache[i] != _cache[i-1]) {
            const FEMCache<FEMCell,FEMVertex>& c=_cache[i];
            FEMInterp I;
            if(c._A->contain(c._B->_pos,I))
                coll.handle(c._bA,c._A,c._bB,c._B,I._coef);
            nrUnique++;
        }
    //INFOV("Nr Cache: %lu, NrUnique: %lu",_cache.size(),nrUnique)
}
void BVHFEMCollision::collideGeom(const FEMGeom& other,FEMCollider& coll,bool surface)
{
    vector<FEMCache<FEMBody,FEMGeomCell> > cache;
    BVHQuery<boost::shared_ptr<FEMBody> > query(*_bvh,dim(),boost::shared_ptr<FEMBody>());
    BVHQuery<boost::shared_ptr<FEMGeomCell> > query2(other.getBVH(),dim(),boost::shared_ptr<FEMGeomCell>());
    query.broadphaseQuery(query2,cache);

    _cacheG.clear();
    std::sort(cache.begin(),cache.end());
    for(sizeType i=0; i<(sizeType)cache.size(); i++)
        if(i == 0 || cache[i] != cache[i-1]) {
            boost::shared_ptr<VBVHFEMBody> bodyA=boost::dynamic_pointer_cast<VBVHFEMBody>(cache[i]._A);
            BVHQuery<sizeType> query(bodyA->_bvh,dim(),bodyA->verbose());
            GeomCallback cb(bodyA,cache[i]._B,_cacheG,dim(),surface);
            query.pointQuery(cb);
        }

    std::sort(_cacheG.begin(),_cacheG.end());
    sizeType nrUnique=0;
    sizeType nrColl=(sizeType)_cacheG.size();
    for(sizeType i=0; i<nrColl; i++)
        if(i == 0 || _cacheG[i] != _cacheG[i-1]) {
            Vec3 n;
            if(_cacheG[i]._B->dist(_cacheG[i]._A->_pos,n))
                coll.handle(_cacheG[i]._bA,_cacheG[i]._A,n);
            nrUnique++;
        }
    //INFOV("Nr Cache: %lu, NrUnique: %lu",_cacheG.size(),nrUnique)
}

//SBVHCollisionBody
class SBVHFEMBody : public VBVHFEMBody
{
public:
    struct SMeshCallback {
        SMeshCallback(boost::shared_ptr<SBVHFEMBody> bA,boost::shared_ptr<SBVHFEMBody> bB,vector<FEMSCache>& cache,sizeType dim,scalar depth)
            :_bA(bA),_bB(bB),_cache(cache),_dim(dim),_depth(depth) {}
        void onCell(const Node<sizeType>& nA,const Node<sizeType>& nB) {
            FEMSCache cache;
            for(sizeType j=0; j<_dim; j++) {
                cache._b[j+1]=_bB;
                cache._v[j+1]=_bB->getVPtr(_bB->_scss[nB._cell][j]);
                cache._n=_bB->_nss[nB._cell];
            }
            cache._depth=_depth;
            for(sizeType i=0; i<_dim; i++) {
                cache._b[0]=_bA;
                cache._v[0]=_bA->getVPtr(_bA->_scss[nA._cell][i]);
                //same element
                for(sizeType j=0; j<_dim; j++)
                    if(cache._v[j+1] == cache._v[0])continue;
                if(nB._bb.distTo(cache._v[0]->_pos) > _depth)continue;
                _cache.push_back(cache);
            }
        }
        boost::shared_ptr<SBVHFEMBody> _bA,_bB;
        vector<FEMSCache>& _cache;
        sizeType _dim;
        scalar _depth;
    };
    SBVHFEMBody():VBVHFEMBody() {}
    virtual FEMBody& operator=(const FEMBody& other) {
        const SBVHFEMBody* otherBVH=static_cast<const SBVHFEMBody*>(&other);
        ASSERT_MSG(otherBVH,"Different Body Type!");
        VBVHFEMBody::operator=(other);
        _scss=otherBVH->_scss;
        return *this;
    }
    bool read(std::istream& is,IOData* dat) {
        VBVHFEMBody::read(is,dat);
        readVector(_scss,is);
        return is.good();
    }
    bool write(std::ostream& os,IOData* dat) const {
        VBVHFEMBody::write(os,dat);
        writeVector(_scss,os);
        return os.good();
    }
    boost::shared_ptr<Serializable> copy() const {
        return boost::shared_ptr<Serializable>(new SBVHFEMBody);
    }
    //build BVH
    virtual void buildBVH() {
        //build faces
        vector<std::pair<Vec3i,Vec2i> > faces;
        getFace(&faces,NULL);
        //build edge
        boost::unordered_map<Vec2i,Vec2i,Hash> edgeMap;
        _scss.resize(faces.size());
        _nss.resize(faces.size());
        _bvh.resize(_scss.size());
        for(sizeType i=0; i<(sizeType)faces.size(); i++) {
            //insert leaf
            _bvh[i]=Node<sizeType>();
            _bvh[i]._nrCell=1;
            _bvh[i]._cell=i;

            //make sure outter normal
            Vec3i& fv=faces[i].first;
            Vec3& n=_nss[i];
            sizeType fo=faces[i].second[1];
            if(dim() == 2) {
                LineSeg seg(_vss[fv[0]]->_pos0,_vss[fv[1]]->_pos0);
                n=seg.normal();
                if(n.dot(_vss[fv[0]]->_pos0-_vss[fo]->_pos0) < 0.0f) {
                    n*=-1.0f;
                    std::swap(fv[0],fv[1]);
                }
                _bvh[i]._bb.setPoints(seg._x,seg._y,seg._y);
            } else {
                Triangle tri(_vss[fv[0]]->_pos0,_vss[fv[1]]->_pos0,_vss[fv[2]]->_pos0);
                n=tri.normal();
                if(n.dot(_vss[fv[0]]->_pos0-_vss[fo]->_pos0) < 0.0f) {
                    n*=-1.0f;
                    std::swap(fv[0],fv[1]);
                }
                _bvh[i]._bb.setPoints(tri._a,tri._b,tri._c);
            }
            _bvh[i]._bb.enlarged(1E-3f,dim());

            //insert surface
            _scss[i].block<3,1>(0,0)=fv;
            _scss[i][3]=faces[i].second[0];
            for(sizeType f=0; f<dim(); f++) {
                Vec2i fid(-1,-1);
                for(sizeType j=0,k=0; j<dim(); j++)
                    if(j!=f)fid[k++]=_scss[i][j];
                std::sort(fid.data(),fid.data()+2);

                if(edgeMap.find(fid) == edgeMap.end())
                    edgeMap[fid]=Vec2i(i,-1);
                else edgeMap[fid][1]=i;
            }
        }
        vector<Face> nodes;
        for(boost::unordered_map<Vec2i,Vec2i,Hash>::const_iterator
                beg=edgeMap.begin(),end=edgeMap.end(); beg!=end; beg++)
            nodes.push_back(Face(Vec3i(beg->first[0],beg->first[1],-1),beg->second));

        //build BVH from bottom up
        sizeType sz=_bvh.size();
        mergeNode(nodes);
        // writeBVHByLevel<sizeType>(_bvh,verbose());
        ASSERT(_bvh.size() == sz*2-1);

        //find optimal depth
        scalar totalV=0.0f;
        for(sizeType i=0; i<(sizeType)_css.size(); i++)
            totalV+=_css[i]->_mass;
        _depth=std::pow(totalV/(scalar)_css.size(),1.0f/(scalar)dim());
        // writeNormal();
    }
    virtual BBox<scalar> updateBVH() {
        sizeType nrSC=(sizeType)_scss.size();
        _nss.resize(_scss.size());
        OMP_PARALLEL_FOR_
        for(sizeType i=0; i<nrSC; i++) {
            _bvh[i]._bb.reset();
            for(sizeType v=0; v<dim(); v++)
                _bvh[i]._bb.setUnion(_vss[_scss[i][v]]->_pos);
            if(dim() == 2)
                _nss[i]=LineSeg(_vss[_scss[i][0]]->_pos,
                                _vss[_scss[i][1]]->_pos).normal();
            else
                _nss[i]=Triangle(_vss[_scss[i][0]]->_pos,
                                 _vss[_scss[i][1]]->_pos,
                                 _vss[_scss[i][2]]->_pos).normal();
        }
        BVHQuery<sizeType>(_bvh,dim(),verbose()).updateBVH();
        return _bvh.back()._bb;
    }
    virtual boost::shared_ptr<FEMCell> getCell(const Node<sizeType>& node) const {
        return _css[_scss[node._cell][3]];
    }
    void writeNormal() {
        updateBVH();
        VTKWriter<scalar> os("Normal","./normal.vtk",true);
        vector<Vec3,Eigen::aligned_allocator<Vec3> > vss;
        for(sizeType i=0; i<(sizeType)_scss.size(); i++) {
            Vec3 pt=Vec3::Zero();
            for(sizeType d=0; d<dim(); d++)
                pt+=_vss[_scss[i][d]]->_pos;
            pt/=(scalar)dim();
            vss.push_back(pt);
            vss.push_back(pt+_nss[i]*_depth);
        }
        os.appendPoints(vss.begin(),vss.end());
        os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,2,0),
                       VTKWriter<scalar>::IteratorIndex<Vec3i>((sizeType)vss.size()/2,2,0),
                       VTKWriter<scalar>::LINE);
    }
    scalar depth() const {
        return _depth;
    }
private:
    vector<Vec4i,Eigen::aligned_allocator<Vec4i> > _scss;
    vector<Vec3,Eigen::aligned_allocator<Vec3> > _nss;
    scalar _depth;
};
boost::shared_ptr<FEMBody> SBVHFEMCollision::createBody() const
{
    return boost::shared_ptr<FEMBody>(new SBVHFEMBody);
}
boost::shared_ptr<FEMCollision> SBVHFEMCollision::copy() const
{
    return boost::shared_ptr<FEMCollision>(new SBVHFEMCollision);
}
void SBVHFEMCollision::collideMesh(FEMCollider& coll,bool surface)
{
    assert(_bvh);
    assert(_mesh);
    vector<FEMCache<FEMBody,FEMBody> > cache;
    BVHQuery<boost::shared_ptr<FEMBody> > query(*_bvh,dim(),boost::shared_ptr<FEMBody>());
    query.broadphaseQuery(query,cache);

    _cacheS.clear();
    std::sort(cache.begin(),cache.end());
    for(sizeType i=0; i<(sizeType)cache.size(); i++)
        if(i == 0 || cache[i] != cache[i-1])
            if(cache[i]._A != cache[i]._B) {
                boost::shared_ptr<SBVHFEMBody> bodyA=boost::dynamic_pointer_cast<SBVHFEMBody>(cache[i]._A);
                boost::shared_ptr<SBVHFEMBody> bodyB=boost::dynamic_pointer_cast<SBVHFEMBody>(cache[i]._B);
				assert(bodyA);
				assert(bodyB);
                BVHQuery<sizeType> query(bodyA->_bvh,dim(),bodyA->verbose());
                BVHQuery<sizeType> query2(bodyB->_bvh,dim(),bodyB->verbose());
                SBVHFEMBody::SMeshCallback cb(bodyA,bodyB,_cacheS,dim(),bodyA->depth());
                query.interBodyDepthQuery(query2,cb,bodyA->depth());
            }

    std::sort(_cacheS.begin(),_cacheS.end());
    sizeType nrColl=(sizeType)_cacheS.size();
    for(sizeType i=0; i<nrColl;) {

        sizeType best=-1;
        scalar bestDist=numeric_limits<scalar>::max(),distTmp;
        Vec3 cp,bestB,bTmp;
		assert(_cacheS[i]._v[0]);
        Vec3 pt=_cacheS[i]._v[0]->_pos;
        //get the list of all collisions with vertex i
        sizeType j=i;
        for(; j<nrColl && _cacheS[j]._v[0] == _cacheS[i]._v[0]; j++)
            if(j == 0 || _cacheS[j] != _cacheS[j-1]) {
                const FEMSCache& c=_cacheS[j];
				assert(c._v[1]);
				assert(c._v[2]);
				assert(c._v[3]);
                if(_mesh->dim() == 2)
                    LineSeg(c._v[1]->_pos,c._v[2]->_pos).calcPointDist(pt,distTmp,cp,bTmp);
                else Triangle(c._v[1]->_pos,c._v[2]->_pos,c._v[3]->_pos).calcPointDist(pt,distTmp,cp,bTmp);
                if(distTmp < bestDist && distTmp < c._depth) {
                    best=j;
                    bestDist=distTmp;
                    bestB=bTmp;
                }
            }

        //check if the vertex is truly inside
        Vec3d coef[4];
        if(best != -1) {
            FEMSCache& c=_cacheS[best];
            if((pt-c._v[1]->_pos).dot(c._n) < 0.0f) {
                coef[0]=c._n.cast<scalarD>();
                coef[1]=coef[0]*-bestB[0];
                coef[2]=coef[0]*-bestB[1];
                coef[3]=coef[0]*-bestB[2];
                coll.handle(c._b,c._v,coef,_mesh->dim()+1);
            }
		  }
          i=j;
        }
    }
