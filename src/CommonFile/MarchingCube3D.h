#ifndef MARCHING_CUBE_3D_H
#define MARCHING_CUBE_3D_H

#include "Config.h"
#include "GridBasic.h"
#include "MarchingCube3DTable.h"
#include "ObjMesh.h"

PRJ_BEGIN

template <typename T>
class MarchingCube3D
{
public:
    typedef typename ScalarUtil<T>::ScalarVec3 PT3;
    MarchingCube3D(ObjMeshTpl<T>& mesh):_mesh(mesh) {}
    virtual ~MarchingCube3D() {}
    void solve(const Grid<T,T>& field,const T& contour,sizeType bdMask=0) {
        BBox<T> bb=field.getBB();
        bb._maxC-=field.getCellSize();
        Grid<Vec3i,T> indexTag;
        indexTag.reset(field.getNrCell()-Vec3i::Ones(),bb,Vec3i::Zero(),field.isCenter());
        Grid<Vec3i,T> vertLookup;
        vertLookup.makeSameGeometry(indexTag);
        solve(field,contour,indexTag,vertLookup,bdMask);
    }
    void solve(const Grid<T,T>& field,const T& contour,
               Grid<Vec3i,T>& indexTag,
               Grid<Vec3i,T>& vertLookup,
               sizeType bdMask=0) {
        ASSERT(indexTag.getBB()._minC == field.getBB()._minC)
        ASSERT((field.getCellSize()-indexTag.getCellSize()).norm() < EPS)
        //ASSERT(!indexTag.isCenter() && !field.isCenter())
        ASSERT(compLE(indexTag.getNrPoint(),field.getNrPoint()-Vec3i::Ones()))
        ASSERT(indexTag.getNrPoint() == vertLookup.getNrPoint())

        Vec3i start(0,0,0);
        Vec3i end=indexTag.getNrPoint();
        Vec3i startVT=start;
        Vec3i endVT=end;
        if(bdMask&1)
            startVT.x()++;
        if(bdMask&2)
            startVT.y()++;
        if(bdMask&4)
            startVT.z()++;
        if(bdMask&8)
            endVT.x()--;
        if(bdMask&16)
            endVT.y()--;
        if(bdMask&32)
            endVT.z()--;
        solve(field,contour,
              indexTag,vertLookup,
              start,end,
              start,end-Vec3i::Ones(),
              startVT,endVT);
    }
    const ObjMeshTpl<T>& getMesh() const {
        return _mesh;
    }
    ObjMeshTpl<T>& getMesh() {
        return _mesh;
    }
    static PT3 genPoint(const PT3& a,const PT3& b,const T& ia,const T& ib) {
        return (b*ia-a*ib)/(ia-ib);
    }
protected:
    void solve(const Grid<T,T>& field,
               const T& contour,
               Grid<Vec3i,T>& indexTag,
               Grid<Vec3i,T>& vertLookup,
               const Vec3i& start,const Vec3i& end,
               const Vec3i& startMC,const Vec3i& endMC,
               const Vec3i& startVT,const Vec3i& endVT) {
        INFO("Start Marching Cube")
        _contour=contour;
        _start=start;
        _end=end;
        _startMC=startMC;
        _endMC=endMC;
        _startVT=startVT;
        _endVT=endVT;
        _nrPoint=indexTag.getNrPoint();
        _stride=field.getStride();
        _strideVert=vertLookup.getStride();
        //offset
        _off[0]=0;
        _off[1]=_stride.z();
        _off[2]=_stride.x()+_stride.z();
        _off[3]=_stride.x();
        _off[4]=_stride.y();
        _off[5]=_stride.z()+_stride.y();
        _off[6]=_stride.x()+_stride.z()+_stride.y();
        _off[7]=_stride.x()+_stride.y();
        //offset vert
        _offVert[0]=0;
        _offVert[1]=_strideVert.z();
        _offVert[2]=_strideVert.x()+_strideVert.z();
        _offVert[3]=_strideVert.x();
        _offVert[4]=_strideVert.y();
        _offVert[5]=_strideVert.z()+_strideVert.y();
        _offVert[6]=_strideVert.x()+_strideVert.z()+_strideVert.y();
        _offVert[7]=_strideVert.x()+_strideVert.y();
        //offset pt
        _offPt[0]=PT3(indexTag.getCellSize().x(),0.0f,0.0f);
        _offPt[1]=PT3(0.0f,indexTag.getCellSize().y(),0.0f);
        _offPt[2]=PT3(0.0f,0.0f,indexTag.getCellSize().z());
        _strideXVert.assign(_nrPoint.x(),0);
        _strideXTri.assign(_nrPoint.x(),0);

        ASSERT(compGE(_start,Vec3i::Zero()))
        ASSERT(compLE(_end,_nrPoint))

        //calculate node type and edge
        INFO("\tAggregating Nodes")
        OMP_PARALLEL_FOR_
        for(sizeType x=_start.x(); x<_end.x(); x++)
            aggregateVerts(field,indexTag,x);
        OMP_PARALLEL_FOR_
        for(sizeType x=_startMC.x(); x<_endMC.x(); x++)
            aggregateTris(field,indexTag,x);

        //scan
        vector<PT3,Eigen::aligned_allocator<PT3> >& vss=_mesh.getV();
        vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& iss=_mesh.getI();
        _totalVert=(sizeType)vss.size();
        _totalTri=(sizeType)iss.size();
        for(sizeType x=0; x<_nrPoint.x(); x++) {
            sizeType tmp=_strideXVert[x];
            _strideXVert[x]=_totalVert;
            _totalVert+=tmp;

            tmp=_strideXTri[x];
            _strideXTri[x]=_totalTri;
            _totalTri+=tmp;
        }
        vss.resize(_totalVert);
        iss.resize(_totalTri);

        //fill vertices
        INFO("\tGenerating Vertices")
        OMP_PARALLEL_FOR_
        for(sizeType x=_startVT.x(); x<_endVT.x(); x++)
            fillVerts(field,indexTag,vertLookup,vss,x,_strideXVert[x]);

        //fill indices
        INFO("\tGenerating Indices")
        OMP_PARALLEL_FOR_
        for(sizeType x=_startMC.x(); x<_endMC.x(); x++)
            fillInds(indexTag,vertLookup,iss,x,_strideXTri[x]);
        INFO("Finish Marching Cube")
    }
    void aggregateVerts(const Grid<T,T>& field,
                        Grid<Vec3i,T>& indexTag,
                        const sizeType& x) {
        //number of vertices
        sizeType& totalVert=_strideXVert[x];

        //calculate how many vert and triangle to generate
        for(sizeType y=_start.y(); y<_end.y(); y++)
            for(sizeType z=_start.z(); z<_end.z(); z++) {
                const sizeType index=field.getIndex(Vec3i(x,y,z));
                Vec3i& val=indexTag.get(Vec3i(x,y,z));
                sizeType& cellType=val.x();
                sizeType& nrVert=val.y();
                sizeType tmp;

                //z  y
                //| /
                //o --x
                cellType=0u;
                cellType+=field.get(index+_off[0]) < _contour ? 1 : 0;
                cellType+=field.get(index+_off[1]) < _contour ? 2 : 0;
                cellType+=field.get(index+_off[2]) < _contour ? 4 : 0;
                cellType+=field.get(index+_off[3]) < _contour ? 8 : 0;
                cellType+=field.get(index+_off[4]) < _contour ? 16 : 0;
                cellType+=field.get(index+_off[5]) < _contour ? 32 : 0;
                cellType+=field.get(index+_off[6]) < _contour ? 64 : 0;
                cellType+=field.get(index+_off[7]) < _contour ? 128 : 0;

                //nr vertex
                bool v0=(cellType&1) != 0;
                bool v1=(cellType&2) != 0;
                bool v3=(cellType&8) != 0;
                bool v4=(cellType&16) != 0;
                nrVert=((v0^v1) ? 1 : 0)+((v0^v3) ? 1 : 0)+((v0^v4) ? 1 : 0);
                tmp=nrVert;
                nrVert=totalVert;
                totalVert+=tmp;
            }
    }
    void aggregateTris(const Grid<T,T>& field,
                       Grid<Vec3i,T>& indexTag,
                       const sizeType& x) {
        //number of vertices
        sizeType& totalTri=_strideXTri[x];

        //calculate how many vert and triangle to generate
        for(sizeType y=_startMC.y(); y<_endMC.y(); y++)
            for(sizeType z=_startMC.z(); z<_endMC.z(); z++) {
                const sizeType index=indexTag.getIndex(Vec3i(x,y,z));
                Vec3i& val=indexTag.get(Vec3i(x,y,z));
                sizeType& cellType=val.x();
                sizeType& nrTri=val.z();
                sizeType tmp;

                //nr triangle
                nrTri=nrTriangles3D[cellType];
                tmp=nrTri;
                nrTri=totalTri;
                totalTri+=tmp;
            }

    }
    void fillVerts(const Grid<T,T>& field,
                   Grid<Vec3i,T>& indexTag,
                   Grid<Vec3i,T>& vertLookup,
                   vector<PT3,Eigen::aligned_allocator<PT3> >& vss,
                   const sizeType& x,const sizeType& vertIndexOff) {
        //generate verts
        for(sizeType y=_startVT.y(); y<_endVT.y(); y++)
            for(sizeType z=_startVT.z(); z<_endVT.z(); z++) {
                const PT3 pt=indexTag.getPt(Vec3i(x,y,z));
                const sizeType index=field.getIndex(Vec3i(x,y,z));
                const Vec3i& val=indexTag.get(Vec3i(x,y,z));
                sizeType cellType=val.x();
                sizeType offVert=val.y()+vertIndexOff;
                Vec3i& vertLU=vertLookup.get(Vec3i(x,y,z));
                vertLU=Vec3i::Constant(-1);

                //nr vertex
                bool v0=(cellType&1) != 0;
                bool v1=(cellType&2) != 0;
                bool v3=(cellType&8) != 0;
                bool v4=(cellType&16) != 0;
                if(v0^v3) {
                    vss[offVert]=
                        genPoint(pt,pt+_offPt[0],
                                 field.get(index)-_contour,
                                 field.get(index+_off[3])-_contour);
                    vertLU.x()=offVert++;
                }
                if(v0^v4) {
                    vss[offVert]=
                        genPoint(pt,pt+_offPt[1],
                                 field.get(index)-_contour,
                                 field.get(index+_off[4])-_contour);
                    vertLU.y()=offVert++;
                }
                if(v0^v1) {
                    vss[offVert]=
                        genPoint(pt,pt+_offPt[2],
                                 field.get(index)-_contour,
                                 field.get(index+_off[1])-_contour);
                    vertLU.z()=offVert++;
                }
            }
    }
    void fillInds(Grid<Vec3i,T>& indexTag,
                  Grid<Vec3i,T>& vertLookup,
                  vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& vss,
                  const sizeType& x,const sizeType& triIndexOff) {
        //generate triangles
        for(sizeType y=_startMC.y(); y<_endMC.y(); y++)
            for(sizeType z=_startMC.z(); z<_endMC.z(); z++) {
                const sizeType index=vertLookup.getIndex(Vec3i(x,y,z));
                const Vec3i& val=indexTag.get(Vec3i(x,y,z));
                sizeType cellType=val.x();
                sizeType offTri=val.z()+triIndexOff;

                sizeType nrTri=nrTriangles3D[cellType];
                for(sizeType idT=0,idOff=cellType*20; idT<nrTri; idT++,idOff+=4) {
                    const sizeType* triId=idTriangles3D+idOff;
                    Vec3i& triRef=vss[offTri++];
                    genTriIndex(triRef.x(),triId[0],index,vertLookup);
                    genTriIndex(triRef.y(),triId[1],index,vertLookup);
                    genTriIndex(triRef.z(),triId[2],index,vertLookup);
                }
            }
    }
    void genTriIndex(sizeType& vertId,const sizeType& edge,const sizeType& index,Grid<Vec3i,T>& vertLookup) const {
        switch(edge) {
        case 0:
            vertId=vertLookup.get(index).z();
            break;
        case 1:
            vertId=vertLookup.get(index+_offVert[1]).x();
            break;
        case 2:
            vertId=vertLookup.get(index+_offVert[3]).z();
            break;
        case 3:
            vertId=vertLookup.get(index).x();
            break;
        case 4:
            vertId=vertLookup.get(index+_offVert[4]).z();
            break;
        case 5:
            vertId=vertLookup.get(index+_offVert[5]).x();
            break;
        case 6:
            vertId=vertLookup.get(index+_offVert[7]).z();
            break;
        case 7:
            vertId=vertLookup.get(index+_offVert[4]).x();
            break;
        case 8:
            vertId=vertLookup.get(index).y();
            break;
        case 9:
            vertId=vertLookup.get(index+_offVert[1]).y();
            break;
        case 10:
            vertId=vertLookup.get(index+_offVert[2]).y();
            break;
        case 11:
            vertId=vertLookup.get(index+_offVert[3]).y();
            break;
        default:
            ASSERT(false);
        }
        ASSERT(vertId != -1);
    }
protected:
    //output
    //_off[0]=;
    //_off[1]=Z;
    //_off[2]=XZ;
    //_off[3]=X;
    //_off[4]=Y;
    //_off[5]=YZ;
    //_off[6]=XYZ;
    //_off[7]=XY;
    sizeType _off[8];
    sizeType _offVert[8];
    PT3 _offPt[3];
    T _contour;
    Vec3i _start;
    Vec3i _end;
    Vec3i _startMC;
    Vec3i _endMC;
    Vec3i _startVT;
    Vec3i _endVT;
    Vec3i _nrPoint;
    Vec3i _stride;
    Vec3i _strideVert;
    std::vector<sizeType> _strideXVert;
    std::vector<sizeType> _strideXTri;
    sizeType _totalVert;
    sizeType _totalTri;
    //result
    ObjMeshTpl<T>& _mesh;
};

PRJ_END

#endif
