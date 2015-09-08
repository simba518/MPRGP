------------------------------------------------------------------------------------------
1. folds

   ./src:              all source codes.
   ./src/MPRGPSolver:  implementation for the traditional and extended MPRGP solver[1], depends on ./src/Utility.
   ./src/Utility: provides a set of auxiliary tools for log, assertion, io, matrix operation, operations on tetrahedron mesh.
   ./src/CommonFile: provides a set of auxiliary tools for ./src/ElasticSolver, implemented by Zherong Pan.
   ./src/ElasticSolver: implements a set of FEM Solver and collision detection, depends on ./scr/CommonFile.
   ./src/CollisionHandle: simulate elastic motions using ./src/ElasticSolver, detect collisions using ./src/ElasticSolver, and resolve collisions using MPRGP solver or ICA solver in ./src/MPRGPSolver.

   ./data: the input data, each model is orgonized by a json file *.ini.

   ./test_case  including a set of test cases for the source code in ./src, as well as a set of data in test_data.


------------------------------------------------------------------------------------------
2. compile
   See CMakeLists.txt. In this version, the CMakeLists.txt only works on Linux. Thus in other systems, you need to build the project by hand.


------------------------------------------------------------------------------------------
3. demo
   see ./script/run.sh
   NOTE 1: it should make sure that the parameter "init_file_dir" in *.ini is correct.
   NOTE 2: for the parameters in *.ini, see void Simulator::init(const string &json_file) in ./src/CollisionHandle/Simulator.cpp


------------------------------------------------------------------------------------------
Authors: Siwang Li, Zherong Pan.
Email: lisiwangcg@gmail.com


------------------------------------------------------------------------------------------
[1] Deformable Objects Collision Handling with Fast Convergence. PG 2015, Siwang Li, Zherong Pan, JinHuang, Xiaogang Jin.
