Utility
=======

A lot of tools for physical simulation.

Installation-------------------------------------------------------------------------
1. dependency

g++
gfortran

boost-dev
boost-system
boost-filesystem

blas/gotoblas/openblas
lapack

superlu2.0
arpack
arpack++

suitesparse
eigen3

opengl
glew
qt4
qglviewer

ipopt
casadi

2. installation

Most of the packages can be installed using synaptic, except the followings:

1) eigen3
   Download the latest version from http://eigen.tuxfamily.org/index.php?title=Main_Page, then follow the instructions in ./INSTALL to install this library.

2) superlu2.0 and arpack++
   arpack++ depends on superlu2.0. However both the superlu2.0 and arpack++ installed from synaptic are broken. I repaired these two packages which can be found in RevisedThirdPartyLibraries on the github. Run ./install to install these two packages.

3) casadi and ipopt
   casadi depends on ipopt. The ipopt library installed from synaptic is broken. Try to compile it from sourse. see https://github.com/casadi/casadi/wiki/InstallationLinux for the installation instructions. 

   Especially, the casadi can be installed by using,
   
   cd ..; mkdir build; cd build
   ../configure --prefix=/usr/local --disable-shared ADD_FFLAGS=-fPIC ADD_CFLAGS=-fPIC ADD_CXXFLAGS=-fPIC
   make; sudo make install