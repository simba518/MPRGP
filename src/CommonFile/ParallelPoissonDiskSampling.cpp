#include "ParallelPoissonDiskSampling.h"
#include "ParticleCD.h"
#include "ObjMesh.h"
#include <set>

USE_PRJ_NAMESPACE
    
ParallelPoissonDiskSampling::ParallelPoissonDiskSampling(sizeType dim)
{
    getDensity()=100.0f;
    getRadius()=0.1f;
    _dim=dim;
    _ratioSz=4;
    _bucketSz=5;
    _k=100;
    _euclid=true;
    buildPhase(4);
}
void ParallelPoissonDiskSampling::sample(ObjMeshTpl<scalar>& mesh)
{
	_euclid=false;
    ASSERT(mesh.getDim() == _dim)
    mesh.smooth();
    generateRawSample(mesh);
    generateHashTable();
    generatePoissonDisk();
}
void ParallelPoissonDiskSampling::sample(
	const std::vector<Vec3,Eigen::aligned_allocator<Vec3> >& vss,
	const std::vector<Vec4i,Eigen::aligned_allocator<Vec4i> >& iss,sizeType* nrP)
{
	_euclid=true;
    generateRawSample(vss,iss,nrP);
    generateHashTable();
    generatePoissonDisk();
}
void ParallelPoissonDiskSampling::sample(const ScalarField& field,bool vol)
{
	_euclid=true;
    if(_dim == 2)
        generateRawSample2D(field,vol);
    else generateRawSample3D(field,vol);
	generateHashTable();
    generatePoissonDisk();
}
void ParallelPoissonDiskSampling::buildPhase(sizeType phaseSz)
{
    _phaseSz=phaseSz;
    BBox<scalar> bb(Vec3::Zero(),Vec3::Ones());
    if(_dim == 2)bb._maxC[2]=0.0f;
    _phase.reset(Vec3i::Constant(_phaseSz),bb,0,true,1);

    sizeType seed=1000;
    Vec3i nrPoint=_phase.getNrPoint();
    sizeType nrPhase=nrPoint.prod();
    _phaseMap.resize(nrPhase);
    std::vector<sizeType> phaseVals(nrPhase);
    for(sizeType i=0;i<nrPhase;i++)
        phaseVals[i]=i;

    for(sizeType x=0;x<nrPoint.x();x++)
    for(sizeType y=0;y<nrPoint.y();y++)
    for(sizeType z=0;z<nrPoint.z();z++){
        sizeType off=std::abs(RandHash::randhashI(seed++))%sizeType(phaseVals.size());
        sizeType val=phaseVals[off];
        _phase.get(Vec3i(x,y,z))=val;
        _phaseMap[val]=PhaseInfo(Vec3i(x,y,z));
        phaseVals.erase(phaseVals.begin()+off);
    }
}
void ParallelPoissonDiskSampling::generateRawSample(ObjMeshTpl<scalar>& mesh)
{
    //generate raw sample
    const std::vector<Vec3,Eigen::aligned_allocator<Vec3> >& vss=mesh.getV();
    const std::vector<Vec3,Eigen::aligned_allocator<Vec3> >& tnss=mesh.getTN();
    const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& iss=mesh.getI();

    const sizeType nrT(iss.size());
    std::vector<scalar> poss(nrT+1);
    poss[0]=0.0f;
    for(sizeType i=1;i<nrT+1;i++)
        poss[i]=poss[i-1]+mesh.getArea(int(i-1));
        
    OMP_PARALLEL_FOR_
    for(sizeType i=0;i<nrT;i++)
        poss[i]/=poss[nrT];
            
    sizeType seed=1000;
    scalar singleArea= _dim == 2 ? _radius*2.0f : _radius*_radius*M_PI;
    _PSet.resize(sizeType(_density*poss[nrT]/singleArea));
            
    //set this to be > 1.0f for safety
    poss[nrT]=2.0f;
    Vec3 bary;
    sizeType tid;
    OMP_PARALLEL_FOR_I(OMP_PRI(bary,tid))
    for(sizeType i=0;i<_PSet.size();i++)
    {
        while(!generateRawPSet(bary,tid,poss,nrT,seed));
        const Vec3i& ids=iss[tid];
        _PSet[i]._pos=vss[ids[0]]*bary[0]+vss[ids[1]]*bary[1];
        if(_dim == 3)
            _PSet[i]._pos+=vss[ids[2]]*bary[2];
        _PSet[i]._normal=tnss[tid];
    }
    _bb=mesh.getBB();
}
void ParallelPoissonDiskSampling::generateRawSample(
	const std::vector<Vec3,Eigen::aligned_allocator<Vec3> >& vss,
	const std::vector<Vec4i,Eigen::aligned_allocator<Vec4i> >& iss,sizeType* nrP)
{
	//generate raw sample
	const sizeType nrT(iss.size());
    std::vector<scalar> poss(nrT+1);
    poss[0]=0.0f;
    for(sizeType i=1;i<nrT+1;i++)
		if(iss[i-1][3] == -1)
        poss[i]=poss[i-1]+TriangleTpl<scalar>
		(vss[iss[i-1][0]],vss[iss[i-1][1]],
		 vss[iss[i-1][2]]).area();
		else
        poss[i]=poss[i-1]+TetrahedronTpl<scalar>
		(vss[iss[i-1][0]],vss[iss[i-1][1]],
		 vss[iss[i-1][2]],vss[iss[i-1][3]]).volume();
	
	if(nrP){
		_radius=poss.back();
		_radius/=(scalar)*nrP;
		if(iss[0][3] == -1)
			_radius=pow((scalar)(_radius/M_PI),(scalar)(1.0f/2.0f));
		else _radius=pow((scalar)(_radius*3.0f/4.0f/M_PI),(scalar)(1.0f/3.0f));
	}

    OMP_PARALLEL_FOR_
    for(sizeType i=0;i<nrT;i++)
        poss[i]/=poss[nrT];
            
    sizeType seed=1000;
    scalar singleArea;
	if(iss[0][3] == -1)
		singleArea=_radius*_radius*M_PI;
	else singleArea=_radius*_radius*_radius*4.0f*M_PI/3.0f;
    _PSet.resize(sizeType(_density*poss[nrT]/singleArea));
            
    //set this to be > 1.0f for safety
	_dim=3;
    poss[nrT]=2.0f;
    Vec4 bary3;
	Vec3 bary2;
    sizeType tid;
	OMP_PARALLEL_FOR_I(OMP_PRI(bary3,bary2,tid))
    for(sizeType i=0;i<_PSet.size();i++)
    {
		if(iss[0][3] == -1){
			while(!generateRawPSet(bary2,tid,poss,nrT,seed));
			const Vec4i& ids=iss[tid];
			_PSet[i]._pos=vss[ids[0]]*bary2[0]+vss[ids[1]]*bary2[1]+
						  vss[ids[2]]*bary2[2];
		}else{
			while(!generateRawPSet(bary3,tid,poss,nrT,seed));
			const Vec4i& ids=iss[tid];
			_PSet[i]._pos=vss[ids[0]]*bary3[0]+vss[ids[1]]*bary3[1]+
						  vss[ids[2]]*bary3[2]+vss[ids[3]]*bary3[3];
		}
		_PSet[i]._normal.setZero();
    }
	_dim=iss[0][3] == -1 ? 2 : 3;
	_bb.reset();
	for(sizeType i=0;i<(sizeType)vss.size();i++)
		_bb.setUnion(vss[i]);
}
void ParallelPoissonDiskSampling::generateRawSample3D(const ScalarField& phi,bool vol)
{
    sizeType seed=1000;
    scalar ratio;
	if(vol)
		ratio=phi.getCellSize().maxCoeff()*
			  phi.getCellSize().maxCoeff()*
			  phi.getCellSize().maxCoeff()/(_radius*_radius*_radius*4.0f*M_PI/3.0f);
	else ratio=phi.getCellSize().maxCoeff()*
			   phi.getCellSize().maxCoeff()/(_radius*_radius*M_PI);
    sizeType nrPCell=(sizeType)std::ceil(scalar(ratio*_density));
    const Vec3i nrPoint=phi.getNrPoint();
    for(sizeType x=0;x<nrPoint.x()-1;x++)
    for(sizeType y=0;y<nrPoint.y()-1;y++)
    for(sizeType z=0;z<nrPoint.z()-1;z++)
    {
        Vec3i id(x,y,z);
        unsigned char inter=
            (phi.get(id+Vec3i(0,0,0)) < 0.0f ? 1 : 0)+
            (phi.get(id+Vec3i(1,0,0)) < 0.0f ? 2 : 0)+
            (phi.get(id+Vec3i(1,1,0)) < 0.0f ? 4 : 0)+
            (phi.get(id+Vec3i(0,1,0)) < 0.0f ? 8 : 0)+
            (phi.get(id+Vec3i(0,0,1)) < 0.0f ? 16 : 0)+
            (phi.get(id+Vec3i(1,0,1)) < 0.0f ? 32 : 0)+
            (phi.get(id+Vec3i(1,1,1)) < 0.0f ? 64 : 0)+
            (phi.get(id+Vec3i(0,1,1)) < 0.0f ? 128 : 0);
        if((vol && inter != 0) || (!vol && inter != 0 && inter != 255))
        {
            for(sizeType i=0;i<nrPCell;i++)
            {
                ParticleN<scalar> p;
                p._pos=phi.getPt(id);
                p._pos[0]+=RandHash::randhash(seed++,0.0f,1.0f)*phi.getCellSize()[0];
                p._pos[1]+=RandHash::randhash(seed++,0.0f,1.0f)*phi.getCellSize()[1];
                p._pos[2]+=RandHash::randhash(seed++,0.0f,1.0f)*phi.getCellSize()[2];
                
				scalar phiVal=phi.sampleSafe3D(p._pos);
                if(!vol || phiVal > 0.0f)
				{
					p._normal=phi.sampleSafe3DGrad(p._pos).normalized();
					p._pos-=p._normal*phiVal;
				}
				_PSet.addParticle(p);
            }
        }
    }
    _bb=phi.getBB();
}
void ParallelPoissonDiskSampling::generateRawSample2D(const ScalarField& phi,bool vol)
{
    sizeType seed=1000;
    scalar ratio;
	if(vol)
		ratio=phi.getCellSize().maxCoeff()*
			  phi.getCellSize().maxCoeff()/(_radius*_radius*M_PI);
	else ratio=phi.getCellSize().maxCoeff()/(_radius*2.0f);
    sizeType nrPCell=(sizeType)std::ceil(scalar(ratio*_density));
    const Vec3i nrPoint=phi.getNrPoint();
    for(sizeType x=0;x<nrPoint.x()-1;x++)
    for(sizeType y=0;y<nrPoint.y()-1;y++)
    {
        Vec3i id(x,y,0);
        unsigned char inter=
            (phi.get(id+Vec3i(0,0,0)) < 0.0f ? 1 : 0)+
            (phi.get(id+Vec3i(1,0,0)) < 0.0f ? 2 : 0)+
            (phi.get(id+Vec3i(1,1,0)) < 0.0f ? 4 : 0)+
            (phi.get(id+Vec3i(0,1,0)) < 0.0f ? 8 : 0);
        if((vol && inter != 0) || (!vol && inter != 0 && inter != 15))
        {
            for(sizeType i=0;i<nrPCell;i++)
            {
                ParticleN<scalar> p;
                p._pos=phi.getPt(id);
                p._pos[0]+=RandHash::randhash(seed++,0.0f,1.0f)*phi.getCellSize()[0];
                p._pos[1]+=RandHash::randhash(seed++,0.0f,1.0f)*phi.getCellSize()[1];
				
				scalar phiVal=phi.sampleSafe2D(p._pos);
                if(!vol || phiVal > 0.0f)
				{
					p._normal=phi.sampleSafe2DGrad(p._pos).normalized();
					p._pos-=p._normal*phiVal;
				}
				_PSet.addParticle(p);
            }
        }
    }
    _bb=phi.getBB();
}
void ParallelPoissonDiskSampling::generateHashTable()
{
    //generate hash table
    {
        if(_dim == 2)
            _cellSz=_radius/sqrt(2.0f);
        else _cellSz=_radius/sqrt(3.0f);
        //very important
        _bb.enlarged(_cellSz*3.0f,_dim);
        Vec3i nrCell=ceil((Vec3)(_bb.getExtent()/_cellSz));
        if(_dim == 2)
            _stride=Vec3i(nrCell[1],1,0);
        else _stride=Vec3i(nrCell[1]*nrCell[2],nrCell[2],1);
    }
    
    //key sorting
    _key.reset(new sizeType[_PSet.size()]);
    _index.reset(new sizeType[_PSet.size()]);
    boost::shared_array<sizeType> keyBK(new sizeType[_PSet.size()]);
    boost::shared_array<sizeType> indexBK(new sizeType[_PSet.size()]);
    OMP_PARALLEL_FOR_
    for(sizeType i=0;i<_PSet.size();i++)
    {
        _key[i]=calcKey(_PSet[i]._pos);
        _index[i]=i;
    }
    CollisionInterface<scalar,ParticleSet>::radixSort(_PSet.size(),_key.get(),keyBK.get(),_index.get(),indexBK.get());
        
    //reorder and sort the phase group
    sizeType nrCell=0;
    sizeType add,phase;
    _rawPSet.resize(_PSet.size());
    for(sizeType i=0;i<sizeType(_phaseMap.size());i++)
        _phaseMap[i]._from=_phaseMap[i]._to=0;
	#ifdef _MSC_VER
	OMP_PARALLEL_FOR_I(OMP_PRI(add,phase),OMP_ADD(nrCell))
    #endif
    for(sizeType i=0;i<_rawPSet.size();i++)
    {
        _rawPSet[i]=_PSet[_index[i]];
        add=(i == 0 || _key[i] != _key[i-1]) ? 1 : 0;
        if(add)
        {
            phase=calcPhase(_rawPSet[i]._pos);
            OMP_CRITICAL_
            _phaseMap[phase]._to++;
        }
        nrCell+=add;
    }
    for(sizeType i=1;i<sizeType(_phaseMap.size());i++)
        _phaseMap[i]._from=_phaseMap[i-1]._from+_phaseMap[i-1]._to;
    for(sizeType i=0;i<sizeType(_phaseMap.size());i++)
        _phaseMap[i]._to=_phaseMap[i]._from;
    _phaseGroup.resize(nrCell);
    OMP_PARALLEL_FOR_I(OMP_PRI(phase))
    for(sizeType i=0;i<_rawPSet.size();i++)
    {
        if(i == 0 || _key[i] != _key[i-1])
        {
            phase=calcPhase(_rawPSet[i]._pos);
            OMP_CRITICAL_
            _phaseGroup[_phaseMap[phase]._to++]=_key[i];
        }
    }

    //generate hash table
    _szTable=nrCell*_ratioSz;
    _hashTable.assign(_szTable*_bucketSz,HashEntry());
    sizeType entry;
    OMP_PARALLEL_FOR_I(OMP_PRI(entry))
    for(sizeType i=0;i<_rawPSet.size();i++)
    {
        if(i == 0){
            entry=insertEntry(_key[i]);
            if(entry >= 0)_hashTable[entry]._from=0;
        }else if(_key[i] != _key[i-1]){
            entry=insertEntry(_key[i]);
            if(entry >= 0)_hashTable[entry]._from=i;
            entry=insertEntry(_key[i-1]);
            if(entry >= 0)_hashTable[entry]._to=i;
        }
        if(i == _rawPSet.size()-1){
            entry=insertEntry(_key[i]);
            if(entry >= 0)_hashTable[entry]._to=i+1;
        }
    }

    //parityCheck();
}
void ParallelPoissonDiskSampling::generatePoissonDisk()
{
    //build neigh list
    std::vector<sizeType> neighs;
    sizeType nrNeigh=0;
    if(_dim == 2){
        for(sizeType x=-2;x<=2;x++)
        for(sizeType y=-2;y<=2;y++){
            Vec3i id(x,y,0);
            if(id != Vec3i::Zero())
            {
                Vec3 minDist=Vec3(std::max<scalar>(scalar(std::abs(x))-1.0f,0.0f),
                                  std::max<scalar>(scalar(std::abs(y))-1.0f,0.0f),0.0f)*_cellSz;
                if(minDist.norm() >= _radius)
                    continue;
                neighs.push_back(id.dot(_stride));
                nrNeigh++;
            }
        }
    }else{
        for(sizeType x=-2;x<=2;x++)
        for(sizeType y=-2;y<=2;y++)
        for(sizeType z=-2;z<=2;z++){
            Vec3i id(x,y,z);
            if(id != Vec3i::Zero())
            {
                Vec3 minDist=Vec3(std::max<scalar>(scalar(std::abs(x))-1.0f,0.0f),
                                  std::max<scalar>(scalar(std::abs(y))-1.0f,0.0f),
                                  std::max<scalar>(scalar(std::abs(z))-1.0f,0.0f))*_cellSz;
                if(minDist.norm() >= _radius)
                    continue;
                neighs.push_back(id.dot(_stride));
                nrNeigh++;
            }
        }
    }

    //assign
    for(sizeType t=0;t<_k;t++)
    {
        for(sizeType p=0;p<sizeType(_phaseMap.size());p++)
        {
            sizeType eid,eido,pid,pido,key;
            OMP_PARALLEL_FOR_I(OMP_PRI(eid,eido,pid,pido,key))
            for(sizeType c=_phaseMap[p]._from;c<_phaseMap[p]._to;c++)
            {
                key=_phaseGroup[c];
                eid=findEntry(key);
                if(eid>=0 && (pid=_hashTable[eid]._from+t) < _hashTable[eid]._to)
                {
                    bool conflict=false;
                    for(sizeType f=0;f<nrNeigh;f++)
                    {
                        eido=findEntry(key+neighs[f]);
                        if(eido>=0 && (pido=_hashTable[eido]._sample) >= 0 && dist(_rawPSet[pid],_rawPSet[pido]) < _radius){
                            conflict=true;
                            break;
                        }
                    }
                    if(!conflict)
                        _hashTable[eid]._sample=pid;
                }
            }
        }
    }

    //assemble
    _PSet.clear();
    for(sizeType i=0;i<sizeType(_hashTable.size());i++)
    {
        if(_hashTable[i]._sample>=0)
            _PSet.addParticle(_rawPSet[_hashTable[i]._sample]);
    }

    //debug
    //for(sizeType i=0;i<_PSet.size();i++)
    //for(sizeType j=i+1;j<_PSet.size();j++)
    //    ASSERT(dist(_PSet[i],_PSet[j]) > _radius)
}
void ParallelPoissonDiskSampling::parityCheck()
{
    //debug particle set sorting
    std::set<sizeType> cell;
    for(sizeType i=0;i<_rawPSet.size();i++)
    {
        ASSERT(_rawPSet[i] == _PSet[_index[i]]);
        cell.insert(calcKey(_rawPSet[i]._pos));
    }
    for(sizeType i=0;i<_rawPSet.size()-1;i++)
        ASSERT(calcKey(_rawPSet[i]._pos) <= calcKey(_rawPSet[i+1]._pos));

    //debug hash table and cell
    sizeType nrCell=0;
    sizeType nrPSet=0;
    for(sizeType i=0;i<sizeType(_hashTable.size());i++)
    {
        HashEntry entry=_hashTable[i];
        if(entry._key >= 0)
        {
            ASSERT(entry._to >= entry._from)
            ASSERT(findEntry(entry._key) == i && entry._sample == -1);
            for(sizeType j=entry._from;j<entry._to;j++){
                ASSERT(calcKey(_rawPSet[j]._pos) == entry._key);
            }
            ASSERT(entry._from == 0 || _key[entry._from-1] != _key[entry._from]);
            ASSERT(entry._to == _rawPSet.size() || _key[entry._to-1] != _key[entry._to]);
            nrCell++;
            nrPSet+=entry._to-entry._from;
        }
    }
    ASSERT(sizeType(cell.size()) == nrCell && _rawPSet.size() == nrPSet);

    //debug phase system
    for(sizeType i=1;i<sizeType(_phaseMap.size());i++)
        ASSERT(_phaseMap[i]._from == _phaseMap[i-1]._to);
    for(sizeType i=0;i<sizeType(_phaseMap.size());i++)
    {
        for(sizeType j=_phaseMap[i]._from;j<_phaseMap[i]._to;j++)
        {
            sizeType entry=findEntry(_phaseGroup[j]);
            ASSERT(entry >= 0);
            ASSERT(_hashTable[entry]._key == _phaseGroup[j]);
        }
    }
}
void findId(sizeType& id,const std::vector<scalar>& poss,const sizeType& nrT,sizeType& seed)
{
	id=0;
    sizeType end=nrT,mid;
    scalar possV=RandHash::randhash(seed++,0.0f,1.0f);
    while(id<end-1)
    {
        mid=(id+end)/2;
        if(poss[mid] < possV)
            id=mid;
        else end=mid;
    }
}
bool ParallelPoissonDiskSampling::generateRawPSet(Vec4& bary,sizeType& id,const std::vector<scalar>& poss,const sizeType& nrT,sizeType& seed) const
{
	findId(id,poss,nrT,seed);
    bary[0]=RandHash::randhash(seed++,0.0f,1.0f);
    bary[1]=RandHash::randhash(seed++,0.0f,1.0f);
    bary[2]=RandHash::randhash(seed++,0.0f,1.0f);
	if(bary[0]+bary[1] > 1.0f){
		bary[0]=1.0f-bary[0];
		bary[1]=1.0f-bary[1];
	}
	if(bary[1]+bary[2] > 1.0f){
		scalar t=bary[2];
		bary[2]=1.0f-bary[0]-bary[1];
		bary[1]=1.0f-t;
	}else if(bary[0]+bary[1]+bary[2] > 1.0f){
		scalar t=bary[2];
		bary[2]=bary[0]+bary[1]+bary[2]-1.0f;
		bary[0]=1.0f-bary[1]-t;
	}
	bary[3]=1.0f-bary[0]-bary[1]-bary[2];
	return (bary[0] >= 0.0f && bary[1] >= 0.0f && bary[2] >= 0.0f && bary[3] >= 0.0f);
}
bool ParallelPoissonDiskSampling::generateRawPSet(Vec3& bary,sizeType& id,const std::vector<scalar>& poss,const sizeType& nrT,sizeType& seed) const
{
    findId(id,poss,nrT,seed);
    bary[0]=RandHash::randhash(seed++,0.0f,1.0f);
    if(_dim == 3){
        bary[1]=sqrt(bary[0])*RandHash::randhash(seed++,0.0f,1.0f);
        bary[0]=1.0f-sqrt(bary[0]);
        bary[2]=1.0f-bary[0]-bary[1];
    }else{
        bary[1]=1.0f-bary[0];
    }
	return (bary[0] >= 0.0f && bary[1] >= 0.0f && bary[2] >= 0.0f);
}
sizeType ParallelPoissonDiskSampling::calcKey(const Vec3& pos) const
{
    Vec3i id=floor((Vec3)((pos-_bb._minC)/_cellSz));
    return id.dot(_stride);
}
sizeType ParallelPoissonDiskSampling::calcPhase(const Vec3& pos) const
{
    Vec3i id=floor((Vec3)((pos-_bb._minC)/_cellSz));
    return _phase.get(Vec3i(id[0]%_phaseSz,id[1]%_phaseSz,id[2]%_phaseSz));
}
sizeType ParallelPoissonDiskSampling::findEntry(const sizeType& key) const
{
    sizeType off=(key%_szTable)*_bucketSz;
    for(sizeType i=0,j=off;i<_bucketSz;i++,j++){
        if(_hashTable[j]._key == -1)
            break;
        else if(_hashTable[j]._key == key)
            return j;
    }
    return -1;
}
sizeType ParallelPoissonDiskSampling::insertEntry(const sizeType& key)
{
    for(sizeType j=(key%_szTable)*_bucketSz,k=j+_bucketSz;j<k;j++)
    {
        OMP_CRITICAL_{
            if(_hashTable[j]._key == -1)
                _hashTable[j]._key=key;
        }
        if(_hashTable[j]._key == key)
            return j;
    }
    return -1;
}
scalar ParallelPoissonDiskSampling::dist(const ParticleN<scalar>& a,const ParticleN<scalar>& b) const
{
    Vec3 dir=b._pos-a._pos;
    scalar de=dir.norm();
    if(_euclid)
        return de;
    
    dir/=de;
    scalar c1=a._normal.dot(dir);
    scalar c2=b._normal.dot(dir);
    if(std::abs(c1-c2) < EPS)
        return de/std::max<scalar>(sqrt(1.0f-c1*c1),0.1f);
    else return de*(asin(c1)-asin(c2))/(c1-c2);
}
