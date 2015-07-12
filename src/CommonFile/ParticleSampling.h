#ifndef PARTICLE_SAMPLING_H
#define PARTICLE_SAMPLING_H

#include "MathBasic.h"
#include "ImplicitFuncUtils.h"

PRJ_BEGIN

class ParticleSampling
{
public:
	static void samplePSet(const ImplicitFunc<scalar>& f,const BBox<scalar>& bb,const Vec3& cellSz,ParticleSetN& solidPSet)
	{
		const Vec3 ext=bb.getExtent();
		const scalar thres=cellSz.maxCoeff()*0.2f;
		const scalar dx=cellSz.maxCoeff()*1.5f;
		const sizeType MAX_ITER=10;
		const Vec3i nrPoint=ceil(Vec3(ext.x()/std::max(cellSz.x(),EPS),
									  ext.y()/std::max(cellSz.y(),EPS),
									  ext.z()/std::max(cellSz.z(),EPS)))+Vec3i::Ones();
		for(sizeType x=0; x<nrPoint.x(); x++)
			for(sizeType y=0; y<nrPoint.y(); y++)
				for(sizeType z=0; z<nrPoint.z(); z++) {
					const Vec3 pt=bb._minC+Vec3((scalar)x*cellSz.x(),(scalar)y*cellSz.y(),(scalar)z*cellSz.z());
					scalar solidPhi=f(pt);
					if(solidPhi < cellSz.maxCoeff() && solidPhi >= -dx*2.0f) {
						ParticleN<scalar> p;
						p._pos=pt;
						//adjust
						sizeType iter=0;
						for(; solidPhi > -thres && iter < MAX_ITER; iter++) {
							Vec3 grad=ImplicitFuncCSGUtils::sampleGrad(f,pt,cellSz).normalized();
							p._pos+=grad*(-thres-solidPhi);
							solidPhi=f(p._pos);
						}
						if(iter == MAX_ITER)
							continue;
						p._vel=ImplicitFuncCSGUtils::sampleGrad(f,pt,cellSz).normalized();
						solidPSet.addParticle(p);
					}
				}
	}
};

PRJ_END

#endif
