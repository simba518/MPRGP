#include <boost/filesystem.hpp>
#include <Eigen/Geometry> 
#include <JsonFilePaser.h>
#include <assertext.h>
#include <FEMGeom.h>
#include <FEMReducedSystem.h>
#include <FEMRigidReducedSystem.h>
#include <FEMWarppedReducedSystem.h>
#include <FEMMesh.h>
#include <ObjMesh.h>
#include <MakeMesh.h>
#include "Simulator.h"
using namespace UTILITY;
USE_PRJ_NAMESPACE

#include <MPRGPSolver.h>
double MATH::ScalarUtil<double>::scalar_max=DBL_MAX;
double MATH::ScalarUtil<double>::scalar_eps=1E-9;

Simulator::Simulator() {

  init_file_name = "./";
  save_results_to = "./";
  time_step = 0.01;
  total_frames = 0;
  invertable = true;
}

void Simulator::init(const string &json_file){

  // open init file
  init_file_name = json_file;
  JsonFilePaser jsonf;
  if ( !jsonf.open(json_file) ){
	return;
  }

  int mesh_group_start_index = 0;
  jsonf.read("invertable", invertable, true);
  
  // set fem solver
  {
	string coll_type = "surface";
	jsonf.read("coll_type", coll_type);
	int coll_type_int = 2; // surface
	if (coll_type == "volume"){
	  coll_type_int = 1;
	}

	string solver_name = "mprgp";
	jsonf.read("fem_solver",solver_name);
	if(solver_name == "penalty"){
	  fem_solver = boost::shared_ptr<FemSolverExt>( new FemSolverExt(coll_type_int) );
	}else if(solver_name == "moseck"){
	  fem_solver = boost::shared_ptr<FemSolverExt>( new MoseckFemSolver(coll_type_int) );
	}else if(solver_name == "ica"){
	  fem_solver = boost::shared_ptr<FemSolverExt>( new ICAFemSolver(coll_type_int) );
	}else{
	  assert(solver_name == "mprgp");
	  fem_solver=boost::shared_ptr<FemSolverExt>(new DecoupledMprgpFemSolver(coll_type_int));
	  // fem_solver = boost::shared_ptr<FemSolverExt>(new MprgpFemSolver(coll_type_int));
	}
  }
  
  // load mesh
  {
	DEBUG_LOG("load mesh");
	vector<string> abq_file;
	if( jsonf.readFilePath("vol_file", abq_file) ){
	  FEMMesh tempt_mesh( fem_solver->getMesh() );
	  for (size_t i = 0; i < abq_file.size(); ++i){
		tempt_mesh.reset(abq_file[i], 0.0f);
		fem_solver->getMesh() += tempt_mesh;
	  }
	}
	vector<vector<double> > trans_rot_scale;
	if(jsonf.read("trans_rot_scale", trans_rot_scale)){
	  for(size_t i = 0; i < trans_rot_scale.size(); i++)
		fem_solver->getMesh().applyTrans(transRotScaleToMat(trans_rot_scale[i]),i,true,true);
	}
  }

  // load mesh group
  {
	DEBUG_LOG("load mesh group");

	mesh_group_start_index = fem_solver->getMesh().nrB();

	vector<string> group_vol_mtl_file;
	vector<double> group_trans_rot_scale_vel;
	vector<double> group_dx_dy_dz_nx_ny_nz;

	if( jsonf.readFilePath("group_vol_mtl_file", group_vol_mtl_file) ){
	  
	  jsonf.read("group_trans_rot_scale_vel", group_trans_rot_scale_vel);
	  jsonf.read("group_dx_dy_dz_nx_ny_nz", group_dx_dy_dz_nx_ny_nz);
	  assert_eq(group_vol_mtl_file.size(), 2);
	  assert_eq(group_trans_rot_scale_vel.size(), 10);
	  assert_eq(group_dx_dy_dz_nx_ny_nz.size(), 6);

	  FEMMesh tempt_mesh( fem_solver->getMesh() );
	  tempt_mesh.reset(group_vol_mtl_file[0], 0.0f);

	  vector<double> trans_rot_scale(7);
	  trans_rot_scale[3] = group_trans_rot_scale_vel[3];
	  trans_rot_scale[4] = group_trans_rot_scale_vel[4];
	  trans_rot_scale[5] = group_trans_rot_scale_vel[5];
	  trans_rot_scale[6] = group_trans_rot_scale_vel[6];

	  Vector3d vel;
	  vel << group_trans_rot_scale_vel[7], group_trans_rot_scale_vel[8], group_trans_rot_scale_vel[9];

	  for (int xi = 0; xi < group_dx_dy_dz_nx_ny_nz[3]; ++xi){
		for (int yi = 0; yi < group_dx_dy_dz_nx_ny_nz[4]; ++yi){
		  for (int zi = 0; zi < group_dx_dy_dz_nx_ny_nz[5]; ++zi){

			// mesh
			const double tx = xi*group_dx_dy_dz_nx_ny_nz[0]+group_trans_rot_scale_vel[0];
			const double ty = yi*group_dx_dy_dz_nx_ny_nz[1]+group_trans_rot_scale_vel[1];
			const double tz = zi*group_dx_dy_dz_nx_ny_nz[2]+group_trans_rot_scale_vel[2];
			fem_solver->getMesh() += tempt_mesh;
			trans_rot_scale[0] = tx;
			trans_rot_scale[1] = ty;
			trans_rot_scale[2] = tz;
			const int body_id = fem_solver->getMesh().nrB()-1; assert_ge(body_id,0);
			fem_solver->getMesh().applyTrans(transRotScaleToMat(trans_rot_scale),body_id,true,true);

			// material
			fem_solver->getMesh().getB(body_id)._system.reset(new FEMSystem(fem_solver->getMesh().getB(body_id)));
			FEMSystem& sys=*(fem_solver->getMesh().getB(body_id)._system);
			sys.clearEnergy();
			sys.readEnergy(group_vol_mtl_file[1], MaterialEnergy::COROTATIONAL, invertable);

			// velocity
			fem_solver->setVel(vel, body_id);
		  }
		}
	  }
	}
  }
  
  // set material
  {
	DEBUG_LOG("set material");
	vector<string> elastic_mtl_file;
	jsonf.readFilePath("elastic_mtl", elastic_mtl_file);
	for (int i = 0; i < fem_solver->getMesh().nrB() && i < mesh_group_start_index; ++i){
	  fem_solver->getMesh().getB(i)._system.reset(new FEMSystem(fem_solver->getMesh().getB(i)));
	  FEMSystem& sys=*(fem_solver->getMesh().getB(i)._system);
	  sys.clearEnergy();
	  if((int)elastic_mtl_file.size() > i){
		sys.readEnergy(elastic_mtl_file[i], MaterialEnergy::COROTATIONAL, invertable);
	  }else{
		sys.addEnergyMaterial(250000.0f, 0.45f, MaterialEnergy::COROTATIONAL, invertable);
	  }
	}

	double ak = 0.01, am = 0.0;
	jsonf.read("alpha_k",ak,0.01);
	jsonf.read("alpha_m",am,0.0);
	assert_ge(ak,0.0f);
	assert_ge(am,0.0f);
	for (int i = 0; i < fem_solver->getMesh().nrB(); ++i){
	  FEMSystem& sys=*(fem_solver->getMesh().getB(i)._system);
	  sys.setDamping(ak, am);
	}

	vector<double> gravity;
	jsonf.read("gravity",gravity);
	assert_eq(gravity.size(),3);
	for (int i = 0; i < fem_solver->getMesh().nrB(); ++i){
	  FEMSystem& sys=*(fem_solver->getMesh().getB(i)._system);
	  sys.addEnergyMass( Vec3(gravity[0], gravity[1], gravity[2]), NULL);
	}

	double coll_k = 1e5;
	jsonf.read("coll_pen", coll_k, 1e5);
	assert_gt(coll_k,0.0f);
	fem_solver->setCollK(coll_k);

	double mu = 0.1f;
	jsonf.read("friction", mu);
	fem_solver->setFriction(mu);
  }

  // set geom
  {
	DEBUG_LOG("set geom");
	boost::shared_ptr<FEMGeom> geom( new FEMGeom(3) );

	// obj mesh
  	vector<string> geom_file;
  	if( jsonf.readFilePath("scene", geom_file) ){

	  vector<vector<double> > scene_trans_scale;
	  jsonf.read("scene_trans_scale", scene_trans_scale);
	  double depth = 0.0;
	  jsonf.read("coll_dectect_depth", depth, 0.0);
	  vector<bool> revert_scene_norm;
	  jsonf.read("revert_scene_norm", revert_scene_norm);

	  for(size_t i = 0; i < geom_file.size(); i++){
		Mat4 R = Mat4::Identity();
		bool revert = false;
		if(scene_trans_scale.size() > i){
		  assert_eq(scene_trans_scale[i].size(), 4);
		  R *= scene_trans_scale[i][3];
		  R(3,3) = 1;
		  R(0,3) = scene_trans_scale[i][0];
		  R(1,3) = scene_trans_scale[i][1];
		  R(2,3) = scene_trans_scale[i][2];
		}
		if(revert_scene_norm.size() > i){
		  revert = revert_scene_norm[i];
		}
		geom->addGeomMesh( R , geom_file[i], depth, revert );
	  }
  	}
	
	// plane
	vector<vector<double> > plane;
	if( jsonf.read("plane_ori_trans_scale", plane) ){
	  Vector4d p;
	  Vector3d trans;
	  for(size_t i = 0; i < plane.size(); i++){
		assert_eq(plane[i].size(),7);
		p << plane[i][0], plane[i][1], plane[i][2], 0;
		trans << plane[i][3], plane[i][4], plane[i][5];
		addPlane(geom, p, trans, plane[i][6]);
	  }
	}

	// sphere
	vector<vector<double> > centre_rad;
	if( jsonf.read("sphere", centre_rad) ){
	  Vector3d ctr;
	  for(size_t i = 0; i < centre_rad.size(); i++){
		assert_eq(centre_rad[i].size(),4);
		ctr << centre_rad[i][0], centre_rad[i][1], centre_rad[i][2];
		geom->addGeomSphere(ctr, centre_rad[i][3]);
	  }
	}
	
	// box
	vector<vector<double> > box_ctr_ext;
	if( jsonf.read("box_ext_trans", box_ctr_ext) ){

	  Vector3d ext;
	  Matrix4d trans = Matrix4d::Identity();
	  for(size_t i = 0; i < box_ctr_ext.size(); i++){
		assert_eq(box_ctr_ext[i].size(),6);
		ext << box_ctr_ext[i][0], box_ctr_ext[i][1], box_ctr_ext[i][2];
		trans(0,3) = box_ctr_ext[i][3];
		trans(1,3) = box_ctr_ext[i][4];
		trans(2,3) = box_ctr_ext[i][5];
		geom->addGeomBox(trans, ext);
	  }
	}

	// cylinder
	vector<vector<double> > cylinder;
	if( jsonf.read("cylinder_rad_y_ori_trans", cylinder) ){
	  Vector4d orient;
	  Vector3d trans;
	  for(size_t i = 0; i < cylinder.size(); i++){
		assert_eq(cylinder[i].size(),8);
		const double rad = cylinder[i][0];
		const double y = cylinder[i][1];
		orient << cylinder[i][2], cylinder[i][3], cylinder[i][4], 0.0;
		trans << cylinder[i][5], cylinder[i][6], cylinder[i][7];
		addCylinder(geom, rad, y, orient, trans, 50, 50);
	  }
	}

	// stair
	vector<double > stair_para;
	if( jsonf.read("stair_ext_trans_yd_zd_n", stair_para) ){
	  assert_eq(stair_para.size(), 9);
	  Vector3d ext;
	  Vector3d trans;
	  ext << stair_para[0], stair_para[1], stair_para[2];
	  trans << stair_para[3], stair_para[4], stair_para[5];
	  const double y_diff = stair_para[6];
	  const double z_diff = stair_para[7];
	  const int num = (int)stair_para[8];
	  addStair(geom, ext, trans, y_diff, z_diff, num);
	}

	// assemble
	geom->assemble();
	fem_solver->_geom = geom;
  }

  // set other parameters
  {
	DEBUG_LOG("set other parameters");
	jsonf.readFilePath("save_to", save_results_to,false);
	fem_solver->setTargetFold(save_results_to);
	jsonf.read("h",time_step,0.01);
	assert_gt(time_step, 0.0f);
	jsonf.read("num_frames", total_frames);
	assert_ge(total_frames, 0);

	int newton_max_iteration = 10, mprgp_max_iteration = 100;
	double newton_tolerance = 1e-4, mprgp_tolerance = 1e-4;

	jsonf.read("newton_max_iteration", newton_max_iteration, 10);
	jsonf.read("newton_tolerance", newton_tolerance, 1e-4);
	jsonf.read("mprgp_max_iteration", mprgp_max_iteration, 100);
	jsonf.read("mprgp_tolerance", mprgp_tolerance, 1e-4);
	
	fem_solver->resetImplicitEuler(newton_tolerance, newton_max_iteration);
	fem_solver->setLinearSolverParameters(mprgp_tolerance, mprgp_max_iteration);

	bool enable_self_con = false;
	jsonf.read("enable_self_con", enable_self_con, false);
	fem_solver->setSelfColl(enable_self_con);

	vector<vector<double> > init_vel;
	if( jsonf.read("init_vel", init_vel) ){
	  Vector3d vel;
	  for (int i = 0; i < fem_solver->getMesh().nrB() && i < (int)init_vel.size() && i < mesh_group_start_index; ++i){
		assert_eq(init_vel[i].size(), 3);
		vel << init_vel[i][0], init_vel[i][1], init_vel[i][2];
		fem_solver->setVel(vel, i);
	  }
	}

	string coll_type = "DCD";
	if( jsonf.read("coll_dect_type",coll_type) && "CCD" == coll_type){
	  fem_solver->setCollDetectionType(FemSolverExt::CCD);
	}else{
	  fem_solver->setCollDetectionType(FemSolverExt::DCD);
	}
  }

  // set planes for collision
  {
	VVec4d planes;
	vector<double> delta_plane;
	vector<vector<double> > plane_constraints;
	if( jsonf.read("plane_con_delta", plane_constraints) ){
	  for (int i = 0; i < plane_constraints.size(); ++i){
		const vector<double> &P = plane_constraints[i];
		assert_eq(P.size(), 5);
		delta_plane.push_back(P[4]);
		Vector4d one_p;
		one_p << P[0], P[1], P[2], P[3];
		planes.push_back(one_p);
	  }
	}
	fem_solver->setCollidePlanes(planes, delta_plane);
  }

  fem_solver->init();

  DEBUG_LOG("finished to load init file.");
}

void Simulator::run(){

  DEBUG_LOG("begin to run simulation.");
  boost::filesystem::create_directory(saveResultsTo());
  boost::filesystem::create_directory(saveResultsTo()+"/binary/");
  boost::filesystem::create_directory(saveResultsTo()+"/abq/");
  boost::filesystem::create_directory(saveResultsTo()+"/QP/");
  boost::filesystem::create_directory(saveResultsTo()+"/collisions/");

  fem_solver->getMesh().writeVTK(saveResultsTo()+"/mesh.vtk");
  fem_solver->_geom->writeVTK(saveResultsTo()+"/scene.vtk");
  for (int i = 0; i < fem_solver->_geom->nrG(); ++i){
	ObjMesh obj_mesh;
    fem_solver->_geom->getG(i).getMesh(obj_mesh);
	ostringstream ossm_scene;
  	ossm_scene << saveResultsTo() << "/scene_" << i << ".obj";
	obj_mesh.write( boost::filesystem::path(ossm_scene.str()) );
  }
  writeAssembledABQ(saveResultsTo()+"/mesh-asembled.abq");

  std::ofstream scene_file(saveResultsTo()+"/binary/scene.b", ios::binary);
  fem_solver->_geom->write(scene_file);

  for (size_t frame = 0; frame < totalFrames(); ++frame){
    
	INFO_LOG("step: " << frame);
	// ostringstream ossm_bin;
  	// ossm_bin << saveResultsTo() << "/binary/frame_" << frame << ".b";
	// std::ofstream mesh_file(ossm_bin.str(), ios::binary);
	// fem_solver->getMesh().write(mesh_file);

	ostringstream ossm_vtk;
  	ossm_vtk << saveResultsTo() << "/frame_" << frame << ".vtk";
  	fem_solver->getMesh().writeVTK( ossm_vtk.str());

	for (int i = 0; i < fem_solver->getMesh().nrB(); ++i){
	  ostringstream ossm_abq;
	  ossm_abq << saveResultsTo() << "/abq/obj_"<<i<< "_frame_" << frame << ".abq";
	  fem_solver->getMesh().getB(i).writeABQ(ossm_abq.str());
	}

  	fem_solver->advance( timeStep() );
  }
}

void Simulator::print()const{
  
  fem_solver->print();
  INFO_LOG("time step: "<< timeStep());
  INFO_LOG("total frames: "<< totalFrames());
  INFO_LOG("save results to: "<< saveResultsTo());
  INFO_LOG("init file:" << init_file_name);
  INFO_LOG("invertable:" << (invertable ? "true":"false"));
}

bool Simulator::writeAssembledABQ(const string filename)const{

  boost::filesystem::ofstream os(filename);

  // write nodes
  os << "*NODE" << std::endl;
  int begin_vindex = 1;
  for (int bi = 0; bi < fem_solver->getMesh().nrB(); ++bi){
	const FEMBody &b = fem_solver->getMesh().getB(bi);
	for (int i = 0; i < (int)b.nrV(); ++i){
	  os<< std::setprecision(12) << (i+begin_vindex) << ", " <<
		b.getV(i)._pos[0] << ", " <<
		b.getV(i)._pos[1] << ", " <<
		b.getV(i)._pos[2] << std::endl;
	}
	begin_vindex += b.nrV();
  }
  
  // write tets
  os << "*ELEMENT, type=C3D4, ELSET=PART1" << std::endl;
  int begin_cindex = 1;
  begin_vindex = 1;
  for (int bi = 0; bi < fem_solver->getMesh().nrB(); ++bi){
	const FEMBody &b = fem_solver->getMesh().getB(bi);
	for (int i = 0; i < (int)b.nrC(); ++i){
	  os<< (i+begin_cindex) << ", " <<
		(b.getC(i)._v[0]->_index+begin_vindex) << ", " <<
		(b.getC(i)._v[1]->_index+begin_vindex) << ", " <<
		(b.getC(i)._v[2]->_index+begin_vindex) << ", " <<
		(b.getC(i)._v[3]->_index+begin_vindex) << std::endl;
	}
	begin_vindex += b.nrV();
	begin_cindex += b.nrC();
  }

  // use same material
  os << "*ELSET, ELSET=EALL, GENERATE" << std::endl;
  os << "1, " << begin_cindex-1 << std::endl;
  return os.good();
}

void Simulator::addStair(boost::shared_ptr<FEMGeom> geom, const Vector3d &ext, Vector3d trans,
						 const double y_diff, const double z_diff, const int num)const{
  
  assert_ge(num, 0);
  assert(geom);

  Matrix4d trans_m = Matrix4d::Identity();
  trans_m(0,3) = trans[0];
  for(int i = 0; i < num; i++){
	trans_m(1,3) = trans[1]+y_diff*(double)i;
	trans_m(2,3) = trans[2]+z_diff*(double)i;
	geom->addGeomBox(trans_m, ext);
  }
}

void Simulator::addPlane(boost::shared_ptr<FEMGeom> geom, const Vector4d &plane,
						 const Vector3d &trans, const double scale)const{
  
  Vector3d ext;
  ext << 20,2,20;
  ext *= scale;
  Matrix4d T = Matrix4d::Identity();
  T.block<3,1>(0,3) = trans;
  geom->addGeomBox(T*orientToTrans(plane), ext);
}

void Simulator::addCylinder(boost::shared_ptr<FEMGeom> geom, 
							const double rad, const double y,
							const Vector4d &orient, const Vector3d &trans,
							const int slice, const int sliceY)const{
  
  ObjMesh mesh;
  MakeMesh::makeCylinder3D(mesh, rad, y, slice, sliceY, true);
  Matrix4d T = Matrix4d::Identity();
  T.block<3,1>(0,3) = trans;
  geom->addGeomMesh(T*orientToTrans(orient), mesh);
}

Matrix4d Simulator::orientToTrans(const Vector4d &orient)const{
  
  const double alpha=-orient[3]/orient.block<3,1>(0,0).squaredNorm();
  const Vector3d p0=orient.block<3,1>(0,0)*alpha;
  Eigen::Quaterniond q;
  q.setFromTwoVectors(Vector3d::Unit(1),orient.block<3,1>(0,0).normalized());

  Matrix4d trans = Matrix4d::Identity();
  trans.block<3,1>(0,3)=p0-orient.block<3,1>(0,0).normalized();
  trans.block<3,3>(0,0)=q.toRotationMatrix();
  return trans;
}

Matrix4d Simulator::transRotScaleToMat(const vector<double> &trans_rot_scale)const{

  assert_eq(trans_rot_scale.size(), 7);

  // scale
  Matrix4d S = Matrix4d::Identity();
  S.block<3,3>(0,0) *= trans_rot_scale[6];

  // rotation
  Matrix4d R = Matrix4d::Identity();
  const double angle_x = trans_rot_scale[3]*M_PI/180.0;
  const double angle_y = trans_rot_scale[4]*M_PI/180.0;
  const double angle_z = trans_rot_scale[5]*M_PI/180.0;
  Matrix3d m;
  m = AngleAxisd(angle_x, Vector3d::UnitX())
	* AngleAxisd(angle_y,  Vector3d::UnitY())
	* AngleAxisd(angle_z, Vector3d::UnitZ());
  R.block<3,3>(0,0) = m;

  // translation
  Matrix4d T = Matrix4d::Identity();
  T.block<3,1>(0,3) << trans_rot_scale[0], trans_rot_scale[1], trans_rot_scale[2];

  // assemble
  return (T*R*S);
}
