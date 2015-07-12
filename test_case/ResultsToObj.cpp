#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>
#include <fstream>
#include <FEMMeshFormat.h>
#include <TetMeshEmbeding.h>
using namespace std;
using namespace boost::filesystem;
using namespace UTILITY;
USE_PRJ_NAMESPACE

void vtkToObj(const string in, const string out){

  
}

VectorXd loadAbqNodes(const string abq_file, const int n){
  
  VectorXd nodes(n);
  boost::filesystem::ifstream is(abq_file);
  int id;
  char comma;
  string line;
  getline(is,line);
  while(is.good()) {
	if(beginsWith(line,"*NODE")) {
	  while(getline(is,line).good() && !line.empty() && line[0] != '*') {
		double x,y,z;
		istringstream iss(line);
		iss >> id >> comma >> x >> comma >> y >> comma >> z;
		nodes[id*3-3] = x;
		nodes[id*3-2] = y;
		nodes[id*3-1] = z;
	  }
	}else getline(is,line);
  }
  return nodes;
}
  
void abqInterpolateObj(pTetMeshEmbeding embed, const VectorXd &rest_shape,
					   const string in, const string out){

  const VectorXd u = loadAbqNodes(in,rest_shape.size())-rest_shape;
  embed->interpolate(u);
  embed->getObjMesh()->write(out);
}

void saveAsObjFiles(const string m_path,const bool interp=false,const string model=""){

  { // abq to obj
	const string abq_fold = m_path+"/abq/";
	const string obj_fold = m_path+"/obj/";
	create_directory(obj_fold);

	pTetMeshEmbeding embed;
	VectorXd rest_shape;
	if(interp){
	  const string rest_tet = model+"/model/mesh.abq";
	  const string rest_obj = model+"/model/mesh.obj";

	  pTetMesh tetMesh = pTetMesh(new TetMesh());
	  pObjmesh objMesh = pObjmesh(new Objmesh());
	  tetMesh->load(rest_tet);
	  objMesh->load(rest_obj);
	  tetMesh->nodes(rest_shape);

	  embed = pTetMeshEmbeding(new TetMeshEmbeding(tetMesh, objMesh));
	  embed->buildInterpWeights();
	}

    directory_iterator end_itr;
	path p (abq_fold.c_str());
    for (directory_iterator itr(p); itr != end_itr; ++itr) {
	  if ( is_regular_file(itr->path()) ) {
		const string file_name = itr->path().filename().string();
		// if (file_name != "obj_0_frame_0.abq")
		//   continue;
		const string abq_file = abq_fold+file_name;
		const string obj_file = obj_fold+file_name.substr(0,file_name.size()-3)+"obj";
		if (boost::filesystem::exists(obj_file)){
		  continue;
		}
		cout << abq_file << endl;
		cout << obj_file << endl;
		if (!interp){
		  FEMMeshFormat::ABQToObj(abq_file, obj_file);
		}else{
		  abqInterpolateObj(embed, rest_shape, abq_file, obj_file);
		}
	  }
	}
  }

}

BOOST_AUTO_TEST_SUITE(ResultsToObj)

BOOST_AUTO_TEST_CASE(ResultsToObj){

  // saveAsObjFiles("./data/dino/tempt_cubes_test",true,"./data/dino/");
  // saveAsObjFiles("./data/dragon/tempt_stair");
  // saveAsObjFiles("./data/bunny/tempt_one2");
  // saveAsObjFiles("./data/longcube/tempt_ica");
  // saveAsObjFiles("./data/longcube/tempt_ica2");
}

BOOST_AUTO_TEST_SUITE_END()
