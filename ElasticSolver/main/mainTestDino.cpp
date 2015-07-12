#include "mainConfig.h"
#ifdef MAIN_TEST_DINO

#include "FEMGeom.h"
#include "FEMMesh.h"
#include "FEMSystem.h"
#include "FEMRigidReducedSystem.h"
#include "ObjMesh.h"
#include <boost/filesystem.hpp>

USE_PRJ_NAMESPACE

//#define USE_REDUCED_SYSTEM
#define FALL_ON_SPHERE
// #define ENABLE_SELF_COLLISION 

const scalar collision_penalty = 1e6f;

void loadMesh(FEMSolver &fem_solver, const string &abq_file)
{
    fem_solver.getMesh().reset(abq_file,0.0f);
}

void setMaterial(FEMSolver &fem_solver,const string& material_path,const string& output)
{
    Vec3 gravity;
    gravity << 0.0f, 0.0f, -9.8f;
#ifdef USE_REDUCED_SYSTEM
    fem_solver.getMesh().getB(0)._system.reset(new FEMRigidReducedSystem(fem_solver.getMesh().getB(0)));
    FEMReducedSystem& sys = (FEMReducedSystem&)*(fem_solver.getMesh().getB(0)._system);
    sys.readEnergy(material_path,MaterialEnergy::LINEAR,true);
    sys.addEnergyMass( gravity, NULL);
    sys.buildU(10,false);
#else
    fem_solver.getMesh().getB(0)._system.reset(new FEMSystem(fem_solver.getMesh().getB(0)));
    FEMReducedSystem& sys = (FEMReducedSystem&)*(fem_solver.getMesh().getB(0)._system);
    sys.readEnergy(material_path,MaterialEnergy::COROTATIONAL,true);
    sys.addEnergyMass( gravity, NULL);
#endif
#ifdef ENABLE_SELF_COLLISION
    fem_solver.setSelfColl(true);
#else
    fem_solver.setSelfColl(false);
#endif
    fem_solver.setCollK(collision_penalty);
    fem_solver.getMesh().getB(0).writeABQ(output);
}

void setGeom(FEMSolver &fem_solver, const string &obj_file_name)
{
    ObjMesh scene_obj;
    std::ifstream obj_file(obj_file_name);
    scene_obj.read(obj_file);
    scene_obj.smooth();
    scene_obj.makeUniform();
    scene_obj.writeVTK("./obj.vtk",true,false,true);

#ifdef FALL_ON_SPHERE
    Mat4 R = Mat4::Identity();
    R.block<3,1>(0,3)=Vec3(-0.5,-0.5,-1.2);
    boost::shared_ptr<FEMGeom> geom( new FEMGeom(3) );
    geom->addGeomMesh( R , scene_obj );

	// add a box as ground
	Vec3 ext;
	ext << 2,2,0.2;
	Mat4 T = Mat4::Identity();
	T.block<3,1>(0,3) << 0,0,-1.4;
	geom->addGeomBox(T, ext);

#else
    Mat4 R = Mat4::Identity();
    R.block<3,3>(0,0).diagonal()*=4.0f;
    R.block<3,1>(0,3)=Vec3(-2.0f,-2.0f,-2.0f);
    boost::shared_ptr<FEMGeom> geom( new FEMGeom(3) );
    geom->addGeomMesh( R , scene_obj, 0.0f, true );
#endif

    geom->assemble();
    fem_solver._geom = geom;
}

void simulateAndSave(FEMSolver &fem_solver, const size_t num_frames, const double time_step, const string &save_to)
{
    boost::filesystem::create_directory(save_to);
    fem_solver.getMesh().writeVTK(save_to+"/mesh.vtk");
    fem_solver._geom->writeVTK(save_to+"/scene.vtk");

    for (size_t frame = 0; frame < num_frames; ++frame) {
        cout << "step: " << frame << endl;
        fem_solver.advance( time_step );
        ostringstream ossm;
        ossm << save_to << "/frame_" << frame << ".vtk";
        fem_solver.getMesh().writeVTK( ossm.str() );
    }
}

int main(int argc, char *argv[])
{
    const string data_root = "./data/dino/";
    const size_t total_frames = 150;
    const double time_step = 0.01f;
    const string abq_file = data_root+"/model/mesh.abq";
    const string scene_obj_file = data_root+"/model/ball.obj";
    const string material_file = data_root+"/model/mesh.elastic";

    FEMSolver fem_solver(3);
    loadMesh(fem_solver, abq_file);
    setMaterial(fem_solver,material_file,abq_file);
    setGeom(fem_solver, scene_obj_file);
    simulateAndSave(fem_solver, total_frames, time_step, data_root+"/tempt/");

    return 0;
}

#endif
