#include "ImplicitFunc.h"

USE_PRJ_NAMESPACE

//ImplicitFuncPlane
ImplicitFuncPlane::ImplicitFuncPlane() {}
ImplicitFuncPlane::ImplicitFuncPlane(const Vec3& x0,const Vec3& n):_p(x0,n) {}
ImplicitFuncPlane::ImplicitFuncPlane(const Vec3& a,const Vec3& b,const Vec3& c):_p(a,b,c) {}
scalar ImplicitFuncPlane::operator()(const Vec3& pos) const
{
    return _p.side(pos);
}
//ImplicitFuncCSG
ImplicitFuncCSG::ImplicitFuncCSG(OP_TYPE op):_a((ImplicitFuncCSG*)NULL),_b((ImplicitFuncCSG*)NULL),_alpha(0.8f),_op(op) {}
ImplicitFuncCSG::ImplicitFuncCSG(const BBox<scalar,3>& bb,OP_TYPE op):_alpha(0.8f),_op(op)
{
    BBox<scalar,2> bb2(bb._minC.block(0,0,2,1),bb._maxC.block(0,0,2,1));
    boost::shared_ptr<ImplicitFuncCSG> axis01(new ImplicitFuncCSG(bb2,op));
    boost::shared_ptr<ImplicitFuncCSG> axis2(new ImplicitFuncCSG(op));
    if(op == INTERSECT) {
        axis2->_a.reset(new ImplicitFuncPlane(Vec3(0.0f,0.0f,bb._minC.z()),Vec3(0.0f,0.0f,-1.0f)));
        axis2->_b.reset(new ImplicitFuncPlane(Vec3(0.0f,0.0f,bb._maxC.z()),Vec3(0.0f,0.0f, 1.0f)));
    } else {
        axis2->_a.reset(new ImplicitFuncPlane(Vec3(0.0f,0.0f,bb._minC.z()),Vec3(0.0f,0.0f, 1.0f)));
        axis2->_b.reset(new ImplicitFuncPlane(Vec3(0.0f,0.0f,bb._maxC.z()),Vec3(0.0f,0.0f,-1.0f)));
    }
    _a=axis01;
    _b=axis2;
}
ImplicitFuncCSG::ImplicitFuncCSG(const BBox<scalar,2>& bb,OP_TYPE op):_alpha(0.8f)
{
    boost::shared_ptr<ImplicitFuncCSG> axis0(new ImplicitFuncCSG(op));
    boost::shared_ptr<ImplicitFuncCSG> axis1(new ImplicitFuncCSG(op));
    if(op == INTERSECT) {
        axis0->_a.reset(new ImplicitFuncPlane(Vec3(bb._minC.x(),0.0f,0.0f),Vec3(-1.0f,0.0f,0.0f)));
        axis0->_b.reset(new ImplicitFuncPlane(Vec3(bb._maxC.x(),0.0f,0.0f),Vec3( 1.0f,0.0f,0.0f)));
        axis1->_a.reset(new ImplicitFuncPlane(Vec3(0.0f,bb._minC.y(),0.0f),Vec3(0.0f,-1.0f,0.0f)));
        axis1->_b.reset(new ImplicitFuncPlane(Vec3(0.0f,bb._maxC.y(),0.0f),Vec3(0.0f, 1.0f,0.0f)));
    } else {
        axis0->_a.reset(new ImplicitFuncPlane(Vec3(bb._minC.x(),0.0f,0.0f),Vec3( 1.0f,0.0f,0.0f)));
        axis0->_b.reset(new ImplicitFuncPlane(Vec3(bb._maxC.x(),0.0f,0.0f),Vec3(-1.0f,0.0f,0.0f)));
        axis1->_a.reset(new ImplicitFuncPlane(Vec3(0.0f,bb._minC.y(),0.0f),Vec3(0.0f, 1.0f,0.0f)));
        axis1->_b.reset(new ImplicitFuncPlane(Vec3(0.0f,bb._maxC.y(),0.0f),Vec3(0.0f,-1.0f,0.0f)));
    }
    _a=axis0;
    _b=axis1;
}
scalar ImplicitFuncCSG::operator()(const Vec3& pos) const
{
    if(!_a && !_b)
        return 1.0f;
    else if(!_b)
        return (*_a)(pos);
    else if(!_a) {
        scalar vb=(*_b)(pos);
        if(_op == SUBTRACT)
            vb*=-1.0f;
        return vb;
    } else {
        scalar va=(*_a)(pos);
        scalar vb=(*_b)(pos);
        if(_op == SUBTRACT)
            vb*=-1.0f;

        scalar sgn=_op == UNION ? -1.0f : 1.0f;
        return (va+vb+sgn*sqrt(std::abs(va*va+vb*vb-2.0f*va*vb*_alpha)))/(1.0f+_alpha);
    }
}
void ImplicitFuncCSG::setAlpha(const scalar& alpha)
{
    _alpha=alpha;
    ImplicitFuncCSG* a=dynamic_cast<ImplicitFuncCSG*>(_a.get());
    ImplicitFuncCSG* b=dynamic_cast<ImplicitFuncCSG*>(_b.get());
    if(a)a->setAlpha(alpha);
    if(b)b->setAlpha(alpha);
}
//ImplicitFuncGridRef
scalar ImplicitFuncGridRef::operator()(const Vec3& pos) const
{
    return _lsRef.sampleSafe(pos);
}
//ImplicitFuncReinit
ImplicitFuncReinit::ImplicitFuncReinit(const Grid<scalar,scalar> &tpl,const ImplicitFunc<scalar>& inner)
{
    _ls.makeSameGeometry(tpl);
    GridOp<scalar,scalar>::copyFromImplictFunc(_ls,inner);
    GridOp<scalar,scalar>::reinitialize(_ls);
}
scalar ImplicitFuncReinit::operator()(const Vec3& pos) const
{
    return _ls.sampleSafe(pos);
}
//ImplicitFuncMesh2D
ImplicitFuncMesh2D::ImplicitFuncMesh2D(const ObjMesh& mesh,const scalar& cellSz,const BBox<scalar>& bb):_cellSz(cellSz)
{
    const vector<Vec3,Eigen::aligned_allocator<Vec3> >& vv=mesh.getV();
    const vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& vi=mesh.getI();
    const BBoxf meshBB=mesh.getBB();
    const scalarF ext=(scalarF)bb.getExtent().maxCoeff();
    _bb=bb;

    _mesh.reset(new MeshSLC());
    //copy verts
    _mesh->_verts.resize(vv.size());
    for(sizeType i=0; i<(sizeType)vv.size(); i++)
        _mesh->_verts[i]=Vec3SLC(vv[i].x(),vv[i].y(),vv[i].z());
    //copy inds
    _mesh->_inds.resize(vi.size());
    for(sizeType i=0; i<(sizeType)vi.size(); i++)
        _mesh->_inds[i]=Vec3i(vi[i].x(),vi[i].y(),vi[i].z());
    //initialize BVH
    _bvh.reset(new AABBvh(_mesh.get()));

    //copy to grid
    Vec3i nrCell=ceil((Vec3)(bb.getExtent()/cellSz));
    BBox<scalar> bbLs(Vec3(bb._minC.x()-ext*0.1f,bb._minC.y()-ext*0.1f,0.0f),
                      Vec3(bb._maxC.x()+ext*0.1f,bb._maxC.y()+ext*0.1f,0.0f));
    _ls.reset(nrCell,bbLs,0.0f);
    for(sizeType x=0; x<nrCell.x(); x++) {
        //INFOV("X: %d",x)
        for(sizeType y=0; y<nrCell.y(); y++) {
            const Vec3 pos=_ls.getPt(Vec3i(x,y,0));
            Vec3SLC a(pos.x(),pos.y(),pos.z());
            a.z()=meshBB._minC.z()-cellSz;
            Vec3SLC b(pos.x(),pos.y(),pos.z());
            b.z()=meshBB._maxC.z()+cellSz;
            _ls.get(Vec3i(x,y,0))=_bvh->intersectLineSeg3D(a,b,0) ? -_cellSz : _cellSz;
        }
    }

    GridOp<scalar,scalar>::reinitialize(_ls);
    GridOp<scalar,scalar>::write2DScalarGridVTK("./levelSet.vtk",_ls);
}
//ImplicitFuncMesh3D
ImplicitFuncMesh3D::ImplicitFuncMesh3D(const ObjMesh& mesh,const scalar& cellSz,const BBox<scalar>& bb):_cellSz(cellSz)
{
    const vector<Vec3,Eigen::aligned_allocator<Vec3> >& vv=mesh.getV();
    const vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& vi=mesh.getI();

    _bb=bb;

    _mesh.reset(new MeshSLC());
    //copy verts
    _mesh->_verts.resize(vv.size());
    for(sizeType i=0; i<(sizeType)vv.size(); i++)
        _mesh->_verts[i]=Vec3SLC(vv[i].x(),vv[i].y(),vv[i].z());
    //copy inds
    _mesh->_inds.resize(vi.size());
    for(sizeType i=0; i<(sizeType)vi.size(); i++)
        _mesh->_inds[i]=Vec3i(vi[i].x(),vi[i].y(),vi[i].z());
    //initialize BVH
    _bvh.reset(new AABBvh(_mesh.get()));
}
scalar ImplicitFuncMesh3D::operator()(const Vec3& pos) const
{
    Vec3SLC a(pos.x(),pos.y(),pos.z());

    //z
    Vec3SLC b(pos.x(),pos.y(),pos.z());
    b.z()=_bb._maxC.z();
    if(!_bvh->intersectLineSeg3D(a,b,0))
        return _cellSz;

    b=Vec3SLC(pos.x(),pos.y(),pos.z());
    b.z()=_bb._minC.z();
    if(!_bvh->intersectLineSeg3D(a,b,0))
        return _cellSz;

    //y
    b=Vec3SLC(pos.x(),pos.y(),pos.z());
    b.y()=_bb._maxC.y();
    if(!_bvh->intersectLineSeg3D(a,b,0))
        return _cellSz;

    b=Vec3SLC(pos.x(),pos.y(),pos.z());
    b.y()=_bb._minC.y();
    if(!_bvh->intersectLineSeg3D(a,b,0))
        return _cellSz;

    //x
    b=Vec3SLC(pos.x(),pos.y(),pos.z());
    b.x()=_bb._maxC.x();
    if(!_bvh->intersectLineSeg3D(a,b,0))
        return _cellSz;

    b=Vec3SLC(pos.x(),pos.y(),pos.z());
    b.x()=_bb._minC.x();
    if(!_bvh->intersectLineSeg3D(a,b,0))
        return _cellSz;

    return -_cellSz;
}
//ImplicitFuncMesh3DAccurate
ImplicitFuncMesh3DAccurate::ImplicitFuncMesh3DAccurate(const ObjMesh& mesh,const scalar& cellSz,const BBox<scalar>& bb):_cellSz(cellSz)
{
    const vector<Vec3,Eigen::aligned_allocator<Vec3> >& vv=mesh.getV();
    const vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& vi=mesh.getI();
    const vector<Vec3,Eigen::aligned_allocator<Vec3> >& tnor=mesh.getTN();
    const vector<Vec3,Eigen::aligned_allocator<Vec3> >& nor=mesh.getN();

    _bb=bb;

    _mesh.reset(new MeshSLC());
    //copy verts
    _mesh->_verts.resize(vv.size());
    for(sizeType i=0; i<(sizeType)vv.size(); i++)
        _mesh->_verts[i]=Vec3SLC(vv[i].x(),vv[i].y(),vv[i].z());
    //copy inds
    _mesh->_inds.resize(vi.size());
    for(sizeType i=0; i<(sizeType)vi.size(); i++)
        _mesh->_inds[i]=Vec3i(vi[i].x(),vi[i].y(),vi[i].z());
    //copy nors
    _mesh->_nors.resize(nor.size());
    for(sizeType i=0; i<(sizeType)nor.size(); i++)
        _mesh->_nors[i]=Vec3SLC(nor[i].x(),nor[i].y(),nor[i].z());
    //copy tnors
    _mesh->_tnors.resize(tnor.size());
    for(sizeType i=0; i<(sizeType)tnor.size(); i++)
        _mesh->_tnors[i]=Vec3SLC(tnor[i].x(),tnor[i].y(),tnor[i].z());
    _mesh->searchNeigh();
    //initialize BVH
    _bvh.reset(new AABBvh(_mesh.get()));
    //initialize octree
    BBox<scalarD> bbD;
    bbD.copy(_bb);
    const Vec3SLC cellSzD(cellSz,cellSz,cellSz);
    sizeType level=0;
    {
        const sizeType maxCellDim=(sizeType)std::ceil(_bb.getExtent().maxCoeff()/cellSz);
        while((((sizeType)1)<<level) < maxCellDim)
            level++;
    }
    INFO("Building Distance Octree");
    _octree.reset(new OctreeDistance(bbD,cellSzD,level,_helper));
    if(!_octree->buildFromAABBvhTopDown(_bvh.get(),false)){
        NOTIFY_MSG("Fail Building OctDistanceTree, Check Your Mesh!")
    }
    INFO("Built");
}
scalar ImplicitFuncMesh3DAccurate::operator()(const Vec3& pos) const
{
    return (scalar)_octree->getPhi(Vec3SLC(pos.x(),pos.y(),pos.z()));
}
//ImplicitFuncRosy
ImplicitFuncRosy::ImplicitFuncRosy(const Vec3& origin,const Vec3& X,const scalar& step,const scalar& coef)
    :_step(step),_coef(coef),_origin(origin),_X(X) {}
scalar ImplicitFuncRosy::operator()(const Vec3& pos) const
{
    Vec3 rel=pos-_origin;
    return _coef*dist(Vec2(rel.dot(_X),(rel-rel.dot(_X)*_X).norm()));
}
scalar ImplicitFuncRosy::dist(const Vec2& p) const
{
    //scalar minX=p.x();
    scalar dist=ScalarUtil<scalar>::scalar_max;
    for(scalar currX=p.x();; currX-=_step) {
        scalar newDist=(Vec2(currX,y(currX))-p).norm();
        if(newDist < dist) {
            dist=newDist;
            //minX=currX;
        } else break;
    }
    for(scalar currX=p.x();; currX+=_step) {
        scalar newDist=(Vec2(currX,y(currX))-p).norm();
        if(newDist < dist) {
            dist=newDist;
            //minX=currX;
        } else break;
    }
    return (y(p.x()) < p.y()) ? dist : -dist;
}
