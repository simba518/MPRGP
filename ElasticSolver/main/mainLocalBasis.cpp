#include "mainConfig.h"
#ifdef MAIN_LOCAL_BASIS

#include "FEMMesh.h"
#include "FEMCollision.h"
#include "FEMLocalBasis.h"
#include "ImplicitFunc.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

class FuncSolid : public ImplicitFunc<scalar>
{
public:
    scalar operator()(const Vec3& pos) const {
        return std::max(std::max(-3.0f-pos.x(),pos.x()-3.0f),
                        std::max(-3.0f-pos.y(),pos.y()-3.0f));
    }
};

int main()
{
    std::string name="CorotationalBeam";
    boost::filesystem::create_directory("./"+name);

    FEMMesh mesh(2,boost::shared_ptr<FEMCollision>(new FEMCollision));
    scalar sz=0.08f;
    BBox<scalar> bb(Vec3(-5.0f,-5.0f,0.0f),Vec3(5.0f,5.0f,0.0f));
    mesh.reset(bb,FuncSolid(),sz,sz*5.0f);
    mesh.writeVTK("./mesh.vtk");

    boost::unordered_map<sizeType,Vec3> constraints;
    FEMLocalBasis basis(mesh.getB(0),constraints);
    basis.setStiffness(0.6f,0.49f);

    FEMBody& body=mesh.getB(0);
    basis.updateMesh(true);
    basis.debugPatchVTK(true,4);
    return 0;
}
#endif
