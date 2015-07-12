#include "FEMLocalBasis.h"
#include "FEMMesh.h"
#include "FEMUtils.h"
#include "FEMSparseReducedBasis.h"
#include "CollisionDetection.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

//helper
FORCE_INLINE sizeType findRoot(vector<sizeType>& loaded,sizeType from)
{
    //find root
    sizeType root=from;
    while(loaded[root] >= 0 && loaded[root] != root)
        root=loaded[root];
    //path compression
    if(root != from)
        loaded[from]=root;
    return root;
}
FEMLocalBasis::FEMLocalBasis(const FEMBody& body,const boost::unordered_map<sizeType,Vec3>& constraints)
    :_body(body),_constraints(constraints)
{
    setStiffness(0.3f,0.45f);
}
void FEMLocalBasis::setStiffness(scalar alpha,scalar poisson)
{
    _a=alpha;
    _b=2.0f*(1.0f-poisson);
}
void FEMLocalBasis::updateMesh(SparseReducedBasis& UL,bool rebuild)
{
    sizeType nrSV=_body.nrSV();
    sizeType nrV=_body.nrV();
    if(rebuild || _ass.size() != nrSV || _fss.size() != nrSV ||
            _conn.rows() != nrV || _conn.cols() != nrV) {
        buildConnectivity();
        buildTruncatedRange();
        buildFrame();
    }
}
void FEMLocalBasis::debugPatchVTK(bool all,sizeType maxRad)
{
    const FixedSparseMatrix<char,Kernel<char> >::ROW& row=_conn.getValue();
    const std::vector<sizeType>& rowOff=_conn.getRowOffset();

    vector<bool> loaded(_body.nrV(),false);
    if(all) {
        for(sizeType c=0; c<(sizeType)loaded.size(); c++)
            loaded[c]=true;
    } else {
        sizeType nrCtr=rand()%_nrSV;
        for(sizeType c=0; c<nrCtr; c++) {
            sizeType cid=rand()%_nrSV;
            sizeType radius=maxRad < 0 ? -maxRad : rand()%(maxRad+1);
            loaded[cid]=true;

            std::deque<sizeType> ss;
            ss.push_back(cid);
            for(sizeType i=0; i<radius; i++) {
                sizeType nrS=ss.size();
                for(sizeType j=0; j<nrS; j++) {
                    sizeType top=ss.front();
                    ss.pop_front();
                    for(sizeType j=rowOff[top]; j<rowOff[top+1] && row[j].second < _nrSV; j++)
                        if(!loaded[row[j].second]) {
                            loaded[row[j].second]=true;
                            ss.push_back(row[j].second);
                        }
                }
            }
        }
    }

    //BlockSparseMatrix UL;
    //generateBasis(loaded,UL,true);
}
sizeType FEMLocalBasis::generateBasis(vector<bool>& loaded,SparseReducedBasis& UL,bool debug) const
{
    return 0;
}
/*sizeType FEMLocalBasis::generateBasis(vector<bool>& loaded,SparseReducedBasis& UL,bool debug) const
{
    //identify each connected patch
    vector<Patch> patches;
    scalar totalA=findPatches(loaded,patches);
    //run K-mean on each patch
    sizeType nrB=0;
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<(sizeType)patches.size(); i++) {
        patches[i]._nrPB=performKMean(patches[i],totalA);
        OMP_CRITICAL_ {
            patches[i]._off=nrB;
            nrB+=patches[i]._nrPB;
        }
    }
    if(debug)
        writePatchVTK(patches,"./patches.vtk");
    //generate basis for each patch
    BASIS_SET basis(nrB);
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<(sizeType)patches.size(); i++) {
        const Patch& p=patches[i];
        for(sizeType v=0; v<(sizeType)p._nodes.size(); v++)
            generateBasis(p._nodes[v][0],p._off+p._nodes[v][1],basis);
        //normalize
        for(sizeType i=0; i<p._nrPB; i++) {
            Vec3d total=Vec3d::Zero();
            for(BASIS::const_iterator beg=basis[p._off+i].begin(),end=basis[p._off+i].end(); beg!=end; beg++)
                for(sizeType d=0; d<_body.dim(); d++)
                    total[d]+=beg->second.col(d).squaredNorm();
            for(sizeType d=0; d<_body.dim(); d++)
                total[d]=sqrt(total[d]);
            for(BASIS::iterator beg=basis[p._off+i].begin(),end=basis[p._off+i].end(); beg!=end; beg++)
                for(sizeType d=0; d<_body.dim(); d++)
                    beg->second.col(d)/=total[d];
        }
    }
    //insert basis into sparse matrix triplets
    for(sizeType i=0; i<(sizeType)basis.size(); i++)
        for(BASIS::const_iterator beg=basis[i].begin(),end=basis[i].end(); beg!=end; beg++)
            addBlock(UL,beg->first*3,i*_body.dim(),beg->second.block(0,0,3,_body.dim()));
    //debug basis
    if(debug)
        writeBasisVTK(basis,"./localBasis/");
    return nrB*_body.dim();
}*/
void FEMLocalBasis::writeVTK(const std::string& path) const
{
    boost::filesystem::create_directory(path);
    writeFrameVTK(path+"/frame.vtk");
    for(sizeType i=0; i<_nrSV; i++) {
        std::ostringstream oss;
        oss << path << "/geodesic" << i << ".vtk";
        writeGeodesicVTK(oss.str(),i);

        std::ostringstream oss2;
        oss2 << path << "/truncatedRange" << i << ".vtk";
        writeTruncatedRangeVTK(oss2.str(),i);
    }
}
//private
scalar FEMLocalBasis::findPatches(vector<bool>& loaded,vector<Patch>& patches) const
{
    //find connected patches
    const FixedSparseMatrix<char,Kernel<char> >::ROW& row=_conn.getValue();
    const std::vector<sizeType>& rowOff=_conn.getRowOffset();

    scalar totalA=0.0f;
    vector<bool> visited(_nrSV,false);
    for(sizeType i=0; i<_nrSV; i++)
        if(!visited[i] && loaded[i]) {
            patches.push_back(Patch());
            Patch& p=patches.back();

            std::deque<sizeType> ss;
            ss.push_back(i);
            visited[i]=true;
            while(!ss.empty()) {
                sizeType top=ss.front();
                ss.pop_front();
                p._area+=_ass[top];
                totalA+=_ass[top];
                p._nodes.push_back(Vec2i(top,0));

                for(sizeType j=rowOff[top]; j<rowOff[top+1] && row[j].second < _nrSV; j++)
                    if(!visited[row[j].second] && loaded[row[j].second]) {
                        visited[row[j].second]=true;
                        ss.push_back(row[j].second);
                    }
            }
        }
    //count surface area and nrPoint
    return totalA;
}
sizeType FEMLocalBasis::performKMean(Patch& p,scalar totalA) const
{
    vector<Vec2i,Eigen::aligned_allocator<Vec2i> >& vss=p._nodes;
    sizeType r=_body._basis ? _body._basis->_U.cols() : 30;
    sizeType nrV=(sizeType)vss.size();
    sizeType nrCtr=std::min<sizeType>((sizeType)std::ceil(p._area*(scalar)r/totalA),nrV);

    vector<Vec4,Eigen::aligned_allocator<Vec4> > ctrs(nrCtr),ctrsNew(nrCtr);
    for(sizeType i=0; i<nrCtr; i++)
        ctrs[i].block<3,1>(0,0)=_body.getV(vss[i][0])._pos0;

    bool changed=true;
    while(changed) {
        changed=false;
        ctrsNew.assign(nrCtr,Vec4::Zero());
        for(sizeType i=0; i<nrV; i++) {
            //decide label
            sizeType lastC=vss[i][1];
            const Vec3& pt=_body.getV(vss[i][0])._pos0;
            scalar dist=numeric_limits<scalar>::max(),distNew;
            for(sizeType c=0; c<nrCtr; c++)
                if((distNew=(pt-ctrs[c].block<3,1>(0,0)).squaredNorm()) < dist) {
                    dist=distNew;
                    vss[i][1]=c;
                }
            if(lastC != vss[i][1])
                changed=true;

            //update ctr
            ctrsNew[vss[i][1]].block<3,1>(0,0)+=pt;
            ctrsNew[vss[i][1]][3]+=1.0f;
        }
        //kickout unused center
        sizeType j=0;
        for(sizeType c=0; c<nrCtr; c++)
            if(ctrsNew[c][3] > 0.0f)
                ctrsNew[j++]=ctrsNew[c];
            else changed=true;
        nrCtr=j;
        //update center
        for(sizeType c=0; c<nrCtr; c++)
            ctrs[c].block<3,1>(0,0)=ctrsNew[c].block<3,1>(0,0)/ctrsNew[c][3];
    }
    return nrCtr;
}
Mat3 FEMLocalBasis::getDisplacement(const Mat3& F,const Vec3& dpt) const
{
    //for the fundamental solution, see:
    //http://en.wikipedia.org/wiki/Linear_elasticity
    Vec3 x=F.transpose()*dpt;
    scalar r=std::pow(x.norm(),_a);
    Mat3 ret=(Mat3::Identity()*(1.0f-1.0f/(2.0f*_b))+1.0f/(2.0f*_b)*(x*x.transpose())/(r*r))/r;
    return (F*ret.transpose())*F.transpose();
}
void FEMLocalBasis::generateBasis(sizeType v,sizeType c,BASIS_SET& basis) const
{
    Vec3 pt0=_body.getV(v)._pos0;
    Mat3 UAvg;

    //iterate through all the neighbors
    const FixedSparseMatrix<char,Kernel<char> >::ROW& row=_conn.getValue();
    const std::vector<sizeType>& rowOff=_conn.getRowOffset();

    sizeType nrV=_body.nrV();
    vector<bool> tag(nrV,false);
    tag[v]=true;

    std::deque<sizeType> ss;
    ss.push_back(v);
    while(!ss.empty()) {
        //find basis
        sizeType top=ss.front();
        ss.pop_front();
        if(_constraints.find(top) != _constraints.end())
            UAvg.setZero();
        else if(top == v) {
            UAvg.setZero();
            scalar WAvg=0.0f;
            for(sizeType j=rowOff[v]; j<rowOff[v+1]; j++)
                if(row[j].second != v) {
                    UAvg+=getDisplacement(_fss[v],_body.getV(row[j].second)._pos0-pt0);
                    WAvg+=1.0f;
                }
            UAvg/=WAvg;
        } else UAvg=getDisplacement(_fss[v],_body.getV(top)._pos0-pt0);

        //truncate basis
        scalar r=(_body.getV(top)._pos0-pt0).norm();
        if(r > _ri) {
            scalar x=(r-_ri)/(_r0-_ri);
            scalar b=4.0f*x*x*x*(1.0f-x)+x*x*x*x;
            UAvg*=1.0f-b;
        }

        //add basis
        BASIS::iterator it=basis[c].find(top);
        if(it == basis[c].end())
            basis[c][top]=UAvg.cast<scalarD>()*_ass[v];
        else it->second+=UAvg.cast<scalarD>()*_ass[v];

        //breadth first search
        for(sizeType j=rowOff[top]; j<rowOff[top+1]; j++)
            if(!tag[row[j].second] && (_body.getV(row[j].second)._pos0-pt0).norm() < _r0) {
                tag[row[j].second]=true;
                ss.push_back(row[j].second);
            }
    }
}
void FEMLocalBasis::buildConnectivity()
{
    sizeType nrV=(sizeType)_body.nrV();
    sizeType nrC=(sizeType)_body.nrC();
    _nrSV=_body.nrSV();
    //build face
    vector<Eigen::Triplet<char,sizeType> > tripsAll;
    for(sizeType i=0; i<nrC; i++) {
        const FEMCell& c=_body.getC(i);
#define ADD_EDGE(I,J)	\
tripsAll.push_back(Eigen::Triplet<char,sizeType>(c._v[I]->_index,c._v[J]->_index,1));	\
tripsAll.push_back(Eigen::Triplet<char,sizeType>(c._v[J]->_index,c._v[I]->_index,1));
        ADD_EDGE(0,1)
        ADD_EDGE(0,2)
        ADD_EDGE(1,2)
        if(_body.dim() == 3) {
            ADD_EDGE(0,3)
            ADD_EDGE(1,3)
            ADD_EDGE(2,3)
        }
    }
    _conn.resize(nrV,nrV);
    _conn.buildFromTripletsDepulicate(tripsAll,0);
    ASSERT(_conn.isSymmetric(0));
}
void FEMLocalBasis::buildTruncatedRange()
{
    sizeType nrV=(sizeType)_body.nrV();
    sizeType nrC=(sizeType)_body.nrC();
    //calculate radius
    sizeType nrS=30;
    if(_body._basis)
        nrS=_body._basis->_U.cols();
    scalar avgMass=0.0f;
    for(sizeType i=0; i<nrC; i++)
        avgMass+=_body.getC(i)._mass;
    avgMass/=(scalar)nrC;
    _r0=std::pow((scalar)(nrS*nrS)*avgMass,1.0f/(scalar)_body.dim());
    _ri=_r0*0.5f;
}
void FEMLocalBasis::buildFrame()
{
    vector<std::pair<Vec3i,Vec2i> > faces;
    _body.getFace(&faces,NULL);

    _ass.assign(_nrSV,0.0f);
    _fss.assign(_nrSV,Mat3::Zero());
    //assemble normal and area
    for(sizeType i=0; i<(sizeType)faces.size(); i++) {
        Vec3i fv=faces[i].first;
        sizeType fo=faces[i].second[1];

        Vec3 n;
        scalar area;
        if(_body.dim() == 2) {
            LineSeg seg(_body.getV(fv[0])._pos0,_body.getV(fv[1])._pos0);
            n=seg.normal();
            if(n.dot(_body.getV(fv[0])._pos0-_body.getV(fo)._pos0) < 0.0f)n*=-1.0f;
            area=seg.length();
        } else {
            Triangle tri(_body.getV(fv[0])._pos0,_body.getV(fv[1])._pos0,_body.getV(fv[2])._pos0);
            n=tri.normal();
            if(n.dot(_body.getV(fv[0])._pos0-_body.getV(fo)._pos0) < 0.0f)n*=-1.0f;
            area=tri.area();
        }
        for(sizeType v=0; v<_body.dim(); v++) {
            _fss[fv[v]].col(1)+=n;
            _ass[fv[v]]+=area/(scalar)_body.dim();
        }
    }
    //build ass and frame
    for(sizeType i=0; i<_nrSV; i++) {
        Mat3& m=_fss[i];
        m.col(1).normalize();
        //find other axis
        if(_body.dim() == 2) {
            m.col(0)=Vec3(m.col(1)[1],-m.col(1)[0],0.0f);
        } else {
            sizeType id;
            m.col(1).cwiseAbs().minCoeff(&id);
            m.col(0)=Vec3::Unit(id);
            m.col(2)=m.col(0).cross(m.col(1)).normalized();
            m.col(0)=m.col(1).cross(m.col(2));
        }
    }
}
//debug
void FEMLocalBasis::writePatchVTK(const vector<Patch>& patches,const std::string& path) const
{
    vector<Vec3,Eigen::aligned_allocator<Vec3> > vss;
    for(sizeType i=0; i<_nrSV; i++)
        vss.push_back(_body.getV(i)._pos0);

    vector<scalar> css(_nrSV,0.0f);
    vector<scalar> lss(_nrSV,0.0f);
    for(sizeType i=0; i<(sizeType)patches.size(); i++)
        for(sizeType pv=0; pv<(sizeType)patches[i]._nodes.size(); pv++) {
            Vec2i n=patches[i]._nodes[pv];
            css[n[0]]=(scalar)(i+1.0f);
            lss[n[0]]=(scalar)(n[1]+1.0f);
        }

    VTKWriter<scalar> os("VTKWriter",path,true);
    os.appendPoints(vss.begin(),vss.end());
    vector<Vec3i,Eigen::aligned_allocator<Vec3i> > iss;
    getSurfaceCss(iss);
    os.appendCells(iss.begin(),iss.end(),_body.dim() == 2 ? VTKWriter<scalar>::LINE : VTKWriter<scalar>::TRIANGLE);
    os.appendCustomPointData("color",css.begin(),css.end());
    os.appendCustomPointData("label",lss.begin(),lss.end());
}
void FEMLocalBasis::writeBasisVTK(const BASIS_SET& basis,const std::string& path) const
{
    sizeType nrC=_body.nrC();
    sizeType nrV=_body.nrV();
    vector<Vec4i,Eigen::aligned_allocator<Vec4i> > iss(nrC);
    for(sizeType i=0; i<nrC; i++)
        for(sizeType j=0; j<_body.dim()+1; j++)
            iss[i][j]=_body.getC(i)(j);
    vector<Vec3,Eigen::aligned_allocator<Vec3> > vss(nrV),vssb;
    for(sizeType v=0; v<nrV; v++)
        vss[v]=_body.getV(v)._pos0;

    boost::filesystem::create_directory(path);
    for(sizeType i=0,k=0; i<(sizeType)basis.size(); i++) {
        for(sizeType j=0; j<_body.dim(); j++,k++) {
            vector<scalar> color(nrV,0.0f);
            vector<scalar> strength(nrV,0.0f);
            vssb=vss;

            for(BASIS::const_iterator beg=basis[i].begin(),end=basis[i].end(); beg!=end; beg++) {
                vssb[beg->first]+=beg->second.col(j).cast<scalar>()*_r0;
                color[beg->first]=1.0f;
                strength[beg->first]=(scalar)(beg->second.col(j).norm());
            }

            std::ostringstream oss;
            oss << path << "/b" << k << ".vtk";
            VTKWriter<scalar> os("Basis",oss.str(),true);
            os.appendPoints(vssb.begin(),vssb.end());
            os.appendCells(iss.begin(),iss.end(),_body.dim() == 2 ? VTKWriter<scalar>::TRIANGLE : VTKWriter<scalar>::TETRA);
            os.appendCustomPointData("range",color.begin(),color.end());
            os.appendCustomPointData("strength",strength.begin(),strength.end());
        }
    }
}
void FEMLocalBasis::writeGeodesicVTK(const std::string& path,sizeType i) const
{
    const FixedSparseMatrix<char,Kernel<char> >::ROW& row=_conn.getValue();
    const std::vector<sizeType>& rowOff=_conn.getRowOffset();

    std::vector<scalar> color(_nrSV,-1.0f);
    std::deque<sizeType> ss;
    ss.push_back(i);
    scalar lv=0.0f;
    while(!ss.empty()) {
        sizeType sz=(sizeType)ss.size();
        for(sizeType i=0; i<sz; i++) {
            sizeType top=ss.front();
            ss.pop_front();
            for(sizeType j=rowOff[top]; j<rowOff[top+1] && row[j].second < _nrSV; j++)
                if(color[row[j].second] == -1.0f) {
                    color[row[j].second]=lv;
                    ss.push_back(row[j].second);
                }
        }
        lv+=1.0f;
    }

    VTKWriter<scalar> os("Geodesic",path,true);
    std::vector<Vec3,Eigen::aligned_allocator<Vec3> > vss(_nrSV);
    for(sizeType i=0; i<_nrSV; i++)
        vss[i]=_body.getV(i)._pos0;
    os.appendPoints(vss.begin(),vss.end());
    vector<Vec3i,Eigen::aligned_allocator<Vec3i> > iss;
    getSurfaceCss(iss);
    os.appendCells(iss.begin(),iss.end(),_body.dim() == 2 ? VTKWriter<scalar>::LINE : VTKWriter<scalar>::TRIANGLE);
    os.appendCustomPointData("Color",color.begin(),color.end());
}
void FEMLocalBasis::writeTruncatedRangeVTK(const std::string& path,sizeType i) const
{
    const FixedSparseMatrix<char,Kernel<char> >::ROW& row=_conn.getValue();
    const std::vector<sizeType>& rowOff=_conn.getRowOffset();

    Vec3 pt0=_body.getV(i)._pos0;
    std::vector<Vec3,Eigen::aligned_allocator<Vec3> > vss;
    std::vector<scalar> color;

    vector<bool> tag(_body.nrV(),false);
    tag[i]=true;

    std::deque<sizeType> ss;
    ss.push_back(i);
    while(!ss.empty()) {
        sizeType top=ss.front();
        ss.pop_front();
        Vec3 pt=_body.getV(top)._pos0;
        vss.push_back(pt);
        color.push_back((pt-pt0).norm());

        for(sizeType j=rowOff[top]; j<rowOff[top+1]; j++)
            if(!tag[row[j].second] && (_body.getV(row[j].second)._pos0-pt0).norm() < _r0) {
                tag[row[j].second]=true;
                ss.push_back(row[j].second);
            }
    }

    VTKWriter<scalar> os("TruncatedRange",path,true);
    os.appendPoints(vss.begin(),vss.end());
    os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,0,0),
                   VTKWriter<scalar>::IteratorIndex<Vec3i>(vss.size(),0,0),VTKWriter<scalar>::POINT);
    os.appendCustomPointData("Color",color.begin(),color.end());
}
void FEMLocalBasis::writeFrameVTK(const std::string& path) const
{
    std::vector<Vec3,Eigen::aligned_allocator<Vec3> > vss;
    std::vector<scalar> color;
    for(sizeType i=0; i<_nrSV; i++) {
        Vec3 pt=_body.getV(i)._pos0;
        for(sizeType d=0; d<_body.dim(); d++) {
            vss.push_back(pt);
            color.push_back((scalar)d);
            vss.push_back(pt+_fss[i].col(d)*_r0*0.1f);
            color.push_back((scalar)d);
        }
    }

    VTKWriter<scalar> os("Frame",path,true);
    os.appendPoints(vss.begin(),vss.end());
    os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,2,0),
                   VTKWriter<scalar>::IteratorIndex<Vec3i>(vss.size()/2,2,0),VTKWriter<scalar>::LINE);
    os.appendCustomPointData("Color",color.begin(),color.end());
}
void FEMLocalBasis::getSurfaceCss(vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& iss) const
{
    for(sizeType i=0; i<_body.nrC(); i++) {
        Vec3i fid;
        const FEMCell& c=_body.getC(i);
        for(sizeType f=0; f<_body.dim()+1; f++) {
            bool surface=true;
            for(sizeType j=0,k=0; j<_body.dim()+1; j++)
                if(j!=f) {
                    fid[k]=c._v[j]->_index;
                    if(fid[k] >= _nrSV) {
                        surface=false;
                        break;
                    }
                    k++;
                }
            if(surface)iss.push_back(fid);
        }
    }
}