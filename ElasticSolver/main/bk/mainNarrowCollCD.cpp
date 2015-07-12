#include "mainConfig.h"
#ifdef MAIN_NARROW_COLLCD_TEST

#include "../FEMMesh.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main()
{
#define DIM 2
#define NR 100
	FEMMesh mesh;
	mesh._dim=DIM;

	boost::filesystem::create_directory("./CollCV");
	for(sizeType i=0;i<NR;i++)
	{
		boost::shared_ptr<FEMMesh::Cell> c(new FEMMesh::Cell);
		for(sizeType d=0;d<DIM+1;d++)
		{
			c->_v[d].reset(new FEMMesh::Vertex);
			c->_v[d]->_pos.setZero();
			c->_v[d]->_pos.block(0,0,DIM,1).setRandom();
		}
		c->_surface=15;
		c->makePositive();
	
		Kernel<scalar>::Vec bary;
		bary.resize(DIM+1);
		bary.setRandom();
		bary.array()=bary.array()*0.5f+1.0f;
		bary/=bary.sum();
		boost::shared_ptr<FEMMesh::Vertex> v(new FEMMesh::Vertex);
	
		v->_pos.setZero();
		for(sizeType d=0;d<DIM+1;d++)
			v->_pos+=(*c)[d]*bary[d];

		ostringstream oss;oss << "./CollCV/f" << i << ".vtk";
		FEMMesh::DebugCollider coll(oss.str());
		mesh.collideCV(c,v,coll);
		
		vector<Vec3,Eigen::aligned_allocator<Vec3> > vss;
		vector<Vec2i,Eigen::aligned_allocator<Vec2i> > iss;
		for(sizeType d=0;d<DIM+1;d++)
			vss.push_back((*c)[d]);
		{
			iss.push_back(Vec2i(0,1));
			iss.push_back(Vec2i(0,2));
			iss.push_back(Vec2i(1,2));
		}
		if(DIM==3)
		{
			iss.push_back(Vec2i(0,3));
			iss.push_back(Vec2i(1,3));
			iss.push_back(Vec2i(2,3));
		}

		ostringstream oss2;oss2 << "./CollCV/tet" << i << ".vtk";
		VTKWriter<scalar> os("Tet",oss2.str(),true);
		os.appendPoints(vss.begin(),vss.end());
		os.appendCells(iss.begin(),iss.end(),VTKWriter<scalar>::LINE);
	}
}
#endif
