#ifndef FEM_MESH_FORMAT_H
#define FEM_MESH_FORMAT_H

#include "MathBasic.h"
#include "ObjMesh.h"

PRJ_BEGIN

class FEMMeshFormat
{
public:
    static void meshToABQ(std::istream& is,std::ostream& os);
    static void meshToABQ(const std::string& is,const std::string& os);
    static void VTKToABQ(std::istream& is,std::ostream& os);
    static void VTKToABQ(const std::string& is,const std::string& os);
    static void ABQToObj(const std::string& is,std::ostream& os);
    static void ABQToObj(const std::string& is,const std::string& os);
    static void segmentObj(scalar deg,ObjMesh& mesh);
    static void genBeam(const Vec3i& nr,scalar cellSz,const std::string& os);
};

PRJ_END

#endif