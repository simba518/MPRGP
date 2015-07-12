#ifndef _SIMULATOR_H_
#define _SIMULATOR_H_

#include <string>
#include "MprgpFemSolver.h"

class Simulator{
	
public:
  Simulator();
  boost::shared_ptr<FemSolverExt> getFemSolver()const{return fem_solver;}
  void init(const string &json_file);
  void run();

  double timeStep()const{
	return time_step;
  }
  size_t totalFrames()const{
	return total_frames;
  }
  const string &saveResultsTo()const{
	return save_results_to;
  }
  void print()const;
  void addStair(boost::shared_ptr<FEMGeom> geom, const Vector3d &ext, Vector3d trans,
				const double height_diff, const double depth_diff, const int num)const;
  void addPlane(boost::shared_ptr<FEMGeom> geom,const Vector4d &plane,
				const Vector3d &trans, const double scale = 1.0f)const;
  void addCylinder(boost::shared_ptr<FEMGeom> geom, 
				   const double rad, const double y,
				   const Vector4d &orient, const Vector3d &trans,
				   const int slice, const int sliceY)const;
  bool writeAssembledABQ(const string filename)const;

protected:
  Matrix4d orientToTrans(const Vector4d &orient)const;
  Matrix4d transRotScaleToMat(const vector<double> &trans_rot_scale)const;

private:
  string init_file_name;
  string save_results_to;
  double time_step;
  size_t total_frames;
  bool invertable;
  boost::shared_ptr<FemSolverExt> fem_solver;
};
  
#endif /*_SIMULATOR_H_*/
