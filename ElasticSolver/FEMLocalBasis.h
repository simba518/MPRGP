#ifndef FEM_LOCAL_BASIS_H
#define FEM_LOCAL_BASIS_H

#include "IO.h"
#include "solvers/MatVec.h"
#include <boost/unordered_map.hpp>

PRJ_BEGIN

class FEMMesh;
struct FEMBody;
struct SparseReducedBasis;
class FEMLocalBasis
{
    typedef vector<Eigen::Triplet<scalarD,sizeType> > TRIPS;
    typedef boost::unordered_map<sizeType,Mat3d> BASIS;
    typedef vector<BASIS> BASIS_SET;
    struct Patch {
        Patch():_area(0.0f) {}
        scalar _area;
        sizeType _off,_nrPB;
        vector<Vec2i,Eigen::aligned_allocator<Vec2i> > _nodes;
    };
public:
    //method
    FEMLocalBasis(const FEMBody& body,const boost::unordered_map<sizeType,Vec3>& constraints);
    void setStiffness(scalar alpha,scalar poisson);
    void updateMesh(SparseReducedBasis& UL,bool rebuild=true);
    void debugPatchVTK(bool all,sizeType maxRad);
    sizeType generateBasis(vector<bool>& loaded,SparseReducedBasis& UL,bool debug=false) const;
    void writeVTK(const std::string& path) const;
private:
    scalar findPatches(vector<bool>& loaded,vector<Patch>& patches) const;
    sizeType performKMean(Patch& p,scalar totalA) const;
    Mat3 getDisplacement(const Mat3& F,const Vec3& dpt) const;
    void generateBasis(sizeType v,sizeType c,BASIS_SET& basis) const;
    void buildConnectivity();
    void buildTruncatedRange();
    void buildFrame();
    //debug
    void writePatchVTK(const vector<Patch>& patches,const std::string& path) const;
    void writeBasisVTK(const BASIS_SET& basis,const std::string& path) const;
    void writeGeodesicVTK(const std::string& path,sizeType i) const;
    void writeTruncatedRangeVTK(const std::string& path,sizeType i) const;
    void writeFrameVTK(const std::string& path) const;
    void getSurfaceCss(vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& iss) const;
    //constraints
    const FEMBody& _body;
    const boost::unordered_map<sizeType,Vec3>& _constraints;
    //data
    FixedSparseMatrix<char,Kernel<char> > _conn;
    vector<Mat3,Eigen::aligned_allocator<Mat3> > _fss;
    vector<scalar> _ass;
    scalar _r0,_ri,_a,_b;
    sizeType _nrSV;
};

PRJ_END

#endif