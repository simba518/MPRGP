#include "mainConfig.h"
#ifdef MAIN_SPARSE_MATRIX

#include "SparseReducedBasis.h"

USE_PRJ_NAMESPACE
	
int main()
{
	BlockSparseMatrix m,mm,mmm;
	
	//build
	m.resize(Vec2i(5,5),Vec2i(10,10));
	Cold blk;
	blk.resize(6);
	blk.setRandom();
	m.addDiagonalBlock(2,7,blk);
	Matd blk2;
	blk2.resize(10,2);
	blk2.setRandom();
	m.addBlock(0,2,blk2);
	m.assemble();
	m.writeVTK("f:/m.vtk","f:/mfrm.vtk");

	Cold x,b;
	x.resize(m.cols());
	b.resize(m.rows());
	x.setRandom();
	m.multiply(x,b);

	Eigen::SparseMatrix<scalarD> sparse;
	sparse.resize((int)m.rows(),(int)m.rows());
	for(int i=0;i<m.rows();i++)
		sparse.coeffRef(i,i)=rand()/(scalarD)RAND_MAX;
	m.leftMultiply(sparse,mm);
	mm.writeVTK("f:/mm.vtk","f:/mmfrm.vtk");

	BlockSparseMatrix::TRIPS trips;
	m.multiplyUTMUAdd(sparse,1.0f,0,0,trips);
	mmm.resize(Vec2i(20,0),Vec2i(20,0));
	mmm.assemble(trips);
	mmm.writeVTK("f:/mmm.vtk","f:/mmmfrm.vtk");
}

#endif