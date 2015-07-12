#ifndef MARCHING_CUBE_2D_H
#define MARCHING_CUBE_2D_H

#include "Config.h"
#include "GridBasic.h"
#include "ObjMesh.h"
#include "MarchingCube3DTable.h"
#include <boost/shared_array.hpp>

PRJ_BEGIN

template <typename T>
class MarchingCube2D
{
public:
    typedef typename ScalarUtil<T>::ScalarVec3 PT3;
    MarchingCube2D(ObjMesh& mesh):_vss(mesh.getV()),_iss(mesh.getI()) {mesh.setDim(2);}
    virtual ~MarchingCube2D() {}
    void solve(const Grid<T,T>& grd,const T& contour=0.0f) {
        _vss.clear();
        _iss.clear();

        ASSERT(grd.getDim() == 2);
        Grid<Vec2i,T> vertIndex;
        vertIndex.makeSameGeometry(grd);
        vertIndex.init(Vec2i::Constant(-1));

        //generate indices
        const Vec3i nrPoint=vertIndex.getNrPoint();
        for(sizeType xx=0; xx<nrPoint.x(); xx++)
            for(sizeType yy=0; yy<nrPoint.y(); yy++) {
                if(xx < nrPoint[0]-1) {
                    if((grd.get(Vec3i(xx,yy,0))-contour)*
                            (grd.get(Vec3i(xx+1,yy,0))-contour) <= 0.0) {
                        vertIndex.get(Vec3i(xx,yy,0))[0]=_vss.size();
                        _vss.push_back(genPoint(grd.getPt(Vec3i(xx,yy,0)),
                                                grd.getPt(Vec3i(xx+1,yy,0)),
                                                grd.get(Vec3i(xx,yy,0))-contour,
                                                grd.get(Vec3i(xx+1,yy,0))-contour));
                    }
                }
                if(yy < nrPoint[1]-1) {
                    if((grd.get(Vec3i(xx,yy,0))-contour)*
                            (grd.get(Vec3i(xx,yy+1,0))-contour) <= 0.0) {
                        vertIndex.get(Vec3i(xx,yy,0))[1]=_vss.size();
                        _vss.push_back(genPoint(grd.getPt(Vec3i(xx,yy,0)),
                                                grd.getPt(Vec3i(xx,yy+1,0)),
                                                grd.get(Vec3i(xx,yy,0))-contour,
                                                grd.get(Vec3i(xx,yy+1,0))-contour));
                    }
                }
            }

        for(int xx=0; xx<nrPoint[0]-1; ++xx)
            for(int yy=0; yy<nrPoint[1]-1; ++yy) {
                char marker=0;
                if(grd.get(Vec3i(xx,yy,0))-contour < 0.0)	 marker+=1;
                if(grd.get(Vec3i(xx+1,yy,0))-contour < 0.0)	 marker+=2;
                if(grd.get(Vec3i(xx+1,yy+1,0))-contour < 0.0)marker+=4;
                if(grd.get(Vec3i(xx,yy+1,0))-contour < 0.0)	 marker+=8;

                sizeType nrTri=nrTriangles2D[(sizeType)marker];
                for(int i=0,j=0; i<nrTri; ++i,j+=2) {
                    Vec3i tmp=Vec3i::Zero();
                    tmp[0]=chooseVert(vertIndex,xx,yy,idTriangles2D[marker*4+j]);
                    tmp[1]=chooseVert(vertIndex,xx,yy,idTriangles2D[marker*4+j+1]);
                    _iss.push_back(tmp);
                }
            }
    }
    static PT3 genPoint(const PT3& a,const PT3& b,const T& ia,const T& ib) {
        PT3 tmp=(b*ia-a*ib)/(ia-ib);
        return PT3(tmp.x(),tmp.y(),0.0f);
    }
    sizeType chooseVert(const Grid<Vec2i,T>& grd,const sizeType& xx,const sizeType& yy,const sizeType& edge) {
        sizeType ret=-1;
        switch(edge) {
        case 0:
            ret=grd.get(Vec3i(xx,yy,0)).y();
            break;
        case 1:
            ret=grd.get(Vec3i(xx+1,yy,0)).y();
            break;
        case 2:
            ret=grd.get(Vec3i(xx,yy,0)).x();
            break;
        case 3:
            ret=grd.get(Vec3i(xx,yy+1,0)).x();
            break;
        }
        ASSERT(ret >= 0);
        return ret;
    }
    std::vector<PT3,Eigen::aligned_allocator<PT3> >& _vss;
    std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& _iss;
};

PRJ_END

#endif
