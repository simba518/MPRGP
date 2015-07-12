#include "mainConfig.h"
#ifdef MAIN_TEST_LOAD_GEOM

#include "FEMGeom.h"
#include "FEMMesh.h"
#include "FEMCollision.h"
#include "FEMCollider.h"
#include "ObjMesh.h"

USE_PRJ_NAMESPACE

void testLoadGeom(const string &geom_file)
{

    cout << "begin of loading: " << geom_file << endl;
    boost::shared_ptr<FEMGeom> geom( new FEMGeom(3) );
    geom->addGeomMesh( Mat4::Identity() , geom_file );
    geom->assemble();
    cout << "end of loading: " << geom_file << endl;
}

int main(int argc, char *argv[])
{

    testLoadGeom("./data/dino/model/ball.obj");
    testLoadGeom("./data/dino/model/plane.obj");
    testLoadGeom("./data/dino/model/cone_head.obj");
    testLoadGeom("./data/dino/model/cube.obj");
    return 0;
}

#endif
