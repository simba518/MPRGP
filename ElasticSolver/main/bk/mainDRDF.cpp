#include "mainConfig.h"
#ifdef MAIN_DRDF_TEST

#include "../RotationUtil.h"
#include <boost/filesystem/operations.hpp>

USE_PRJ_NAMESPACE

int main()
{
	for(sizeType test=0;test<100;test++)
	{
		Mat3 F=Mat3::Random();
		Eigen::Matrix<scalarD,9,9> dRdF;

		Eigen::JacobiSVD<Mat3> svdF(F,Eigen::ComputeFullU|Eigen::ComputeFullV);
		derivRDF(dRdF,F,svdF.matrixU(),svdF.singularValues(),svdF.matrixV());
		Mat3 R=svdF.matrixU()*svdF.matrixV().transpose();

#define DELTA 1E-7f
		for(sizeType r=0;r<3;r++)
		for(sizeType c=0;c<3;c++)
		{
			scalar tmp=F(r,c);
			F(r,c)+=DELTA;
		
			Eigen::JacobiSVD<Mat3> svdFTmp(F,Eigen::ComputeFullU|Eigen::ComputeFullV);
			Mat3 R2=svdFTmp.matrixU()*svdFTmp.matrixV().transpose();
			R2=(R2-R)/DELTA;
			INFOV("%f %f",dRdF.col(r+c*3).norm(),R2.norm())
			R2-=Eigen::Map<Mat3d>(dRdF.col(r+c*3).data()).cast<scalar>();
			INFOV("Err: %f",R2.norm())

			F(r,c)=tmp;
		}
#undef DELTA
	}
}
#endif
