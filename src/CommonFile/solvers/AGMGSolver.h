#ifndef AGMG_SOLVER_H
#define AGMG_SOLVER_H

#include "solvers/AMGSolver.h"

PRJ_BEGIN

//Aggregation based multigrid-K-Cycle solver, highly optimized version
template <typename T>
class AGMGSolver : public AMGSolver<T>
{
public:
    struct AggregationMatrix {
        std::vector<sizeType> _parent;
        sizeType _nrRow;
    };
    struct AGMGLevel : public AMGLevel {
        virtual void getMemoryConsumption(T& MBA,T& MBXRVD) const {
            _A.getMemoryConsumption(MBA);
            MBXRVD =(T)sizeof(T)*(T)(_x.size()+_r.size()+_v.size()+_d.size());
            MBXRVD+=(T)sizeof(sizeType)*(T)_parent._parent.size();
            MBXRVD/=(1024.0f*1024.0f);
        }
        AggregationMatrix _parent;
        Vec _v;
        Vec _d;
    };
    AGMGSolver(bool FCG=true) {
        _nrKCycle=1;
        _FCG=FCG;
        _checkDD=true;
        _efficient=true;
        _vCycle=false;
		_reuseAggr=false;
    }
	void setReuseAggr(bool reuse){_reuseAggr=reuse;}
    virtual void setMatrix(const FixedSparseMatrix<T,KERNEL_TYPE>& matrix,bool syncPrecon) {
        FixedSparseMatrix<T,KERNEL_TYPE> currentMatrix=matrix;
        boost::shared_ptr<AGMGLevel> currentLevel;

		if(_reuseAggr)
		{
			ASSERT(_A == &matrix);
			for(sizeType i=1;i<(sizeType)_levels.size();i++)
			{
				AGMGLevel& lv=(AGMGLevel&)*(_levels[i]);
				buildLevelMatrix(lv._parent,currentMatrix);
				lv._A=currentMatrix;
			}
		}
		else
		{
			//build base level
			_nrLv=0;
			currentLevel.reset((AGMGLevel*)NULL);
			_A=&matrix;

			if((sizeType)_levels.size() <= _nrLv)
				_levels.resize(_nrLv+1);
			_levels[_nrLv]=currentLevel;
			_nrLv++;

			//coarsen
			sizeType lastNrC=currentMatrix.rows();
			while(true) {
				currentLevel.reset(getNewLevel());
				aggregate(currentMatrix,currentLevel->_parent);
				sizeType nrC=currentLevel->_parent._nrRow;
				if(nrC <= 1 || nrC >= lastNrC*0.75f)
					break;
				else {
					buildLevelMatrix(currentLevel->_parent,currentMatrix);
					currentLevel->_A=currentMatrix;

					if((sizeType)_levels.size() <= _nrLv)
						_levels.resize(_nrLv+1);
					_levels[_nrLv]=currentLevel;
					_nrLv++;
				}
				lastNrC=nrC;
			}
		}
        initMem(true);
		debugFillin();
    }
    virtual void setMatrix(const FixedSparseMatrix<T,KERNEL_TYPE>& A0,std::vector<FixedSparseMatrix<T,KERNEL_TYPE> >& ITs,bool smooth) {
        FixedSparseMatrix<T,KERNEL_TYPE> lastA=A0;
        _nrLv=(sizeType)ITs.size()+1;
        if((sizeType)_levels.size() < _nrLv)
            _levels.resize(_nrLv);

        //base level
        _levels[0].reset((AMGLevel*)NULL);
        _A=&A0;

        //build coarser levels
        for(sizeType i=0; i<(sizeType)ITs.size(); i++) {
            _levels[i+1].reset(getNewLevel());
            //build the coarse level matrix
            buildParent(ITs[i],((AGMGLevel&)*(_levels[i+1]))._parent);
            buildLevelMatrix(((AGMGLevel&)*(_levels[i+1]))._parent,lastA);
            _levels[i+1]->_A=lastA;
        }
        initMem(true);
		debugFillin();
    }
protected:
    //aggregate and level build
    void aggregate(const FixedSparseMatrix<T,KERNEL_TYPE>& sm,AggregationMatrix& parent) {
        FixedSparseMatrix<T,KERNEL_TYPE> smTmp=sm;
        aggregatePairwise(sm,parent,_checkDD);
        if(parent._nrRow <= 1)
			return;
		buildLevelMatrix(parent,smTmp);

        AggregationMatrix parent2;
        aggregatePairwise(smTmp,parent2,false);
        sizeType nrC=(sizeType)parent._parent.size();
        for(sizeType i=0; i<nrC; i++) {
            sizeType& p=parent._parent[i];
            if(p >= 0 && parent2._parent[p] >= 0)
                p=parent2._parent[p];
            else p=-1;
        }
        parent._nrRow=parent2._nrRow;
    }
    void aggregatePairwise(const FixedSparseMatrix<T,KERNEL_TYPE>& sm,AggregationMatrix& parent,bool checkDD) const {
        std::vector<sizeType> heap;
        std::vector<sizeType> heapOff(sm.rows());
        std::vector<sizeType> lambda(sm.rows(),0);
        std::vector<std::set<sizeType> > S(sm.rows());
        std::vector<T> maxCouple(sm.rows());
        std::vector<unsigned char> valid(sm.rows(),true);
        parent._parent.resize(sm.rows());

        //coarsen phase I
        for(sizeType i=0; i<sm.rows(); i++) {
            T maxVal=0.0f;
            T total=0.0f;
            for(ConstSMIterator<T> beg=sm.begin(i),end=sm.end(i); beg!=end; ++beg) {
                if(beg.col() != i) {
                    maxVal=std::max(maxVal,-(*beg));
                    total+=std::abs((*beg));
                }
            }
            maxCouple[i]=maxVal;
            if(checkDD && sm(i,i) > 5.0f*total)
                valid[i]=false;
            else if(sm(i,i) < ScalarUtil<T>::scalar_eps)
                valid[i]=false;

            for(ConstSMIterator<T> beg=sm.begin(i),end=sm.end(i); beg!=end; ++beg) {
                if(beg.col() != i && -(*beg) >= maxCouple[i]*_alpha) {
                    lambda[beg.col()]++;
                    S[i].insert(beg.col());
                }
            }
        }
        for(sizeType i=0; i<sm.rows(); i++) {
            if(valid[i])
                pushHeapDef(lambda,heapOff,heap,i);
        }
        sizeType idR=0;
        while(!heap.empty()) {
            sizeType err;
            sizeType i=popHeapDef(lambda,heapOff,heap,err);
            sizeType j=-1;
            if(!valid[i])
                continue;

            //find the other element
            T minCoef=ScalarUtil<T>::scalar_max;
            for(std::set<sizeType>::const_iterator beg=S[i].begin(),end=S[i].end(); beg!=end; ++beg) {
                T coef=sm(i,*beg);
                if(valid[*beg] && coef < minCoef) {
                    j=*beg;
                    minCoef=coef;
                }
            }
            //update insert pair
            valid[i]=0;
            parent._parent[i]=idR;
            if(j!=-1) {
                valid[j]=0;
                parent._parent[j]=idR;
            }
            idR++;
            //update lambda set
            for(std::set<sizeType>::const_iterator beg=S[i].begin(),end=S[i].end(); beg!=end; ++beg) {
                if(valid[*beg]) {
                    lambda[*beg]--;
                    updateHeapDef(lambda,heapOff,heap,*beg);
                }
            }
            if(j!=-1)
                for(std::set<sizeType>::const_iterator beg=S[j].begin(),end=S[j].end(); beg!=end; ++beg) {
                    if(valid[*beg]) {
                        lambda[*beg]--;
                        updateHeapDef(lambda,heapOff,heap,*beg);
                    }
                }
        }
        parent._nrRow=idR;
    }
    void buildParent(const FixedSparseMatrix<T,Kernel<T> >& IT,AggregationMatrix& parent) const {
        const std::vector<sizeType>& rowStart=IT.getRowOffset();
        const FixedSparseMatrix<T,KERNEL_TYPE>::ROW& value=IT.getValue();

        parent._nrRow=IT.rows();
        parent._parent.assign(IT.cols(),-1);
        for(sizeType i=0; i<IT.rows(); i++) {
            for(sizeType j=rowStart[i]; j<rowStart[i+1]; ++j)
                if(abs(value[j].first-1.0f) < 1E-5f)
                    parent._parent[value[j].second]=i;
        }
    }
    void buildLevelMatrix(const AggregationMatrix& parent,FixedSparseMatrix<T,Kernel<T> >& M) const {
        //compute IT*A*I
        const std::vector<sizeType>& rowStart=M.getRowOffset();
        const FixedSparseMatrix<T,KERNEL_TYPE>::ROW& value=M.getValue();

        std::vector<Eigen::Triplet<T,sizeType> > trips;
        trips.reserve(value.size());
        for(sizeType i=0; i<M.rows(); i++) {
            for(sizeType j=rowStart[i]; j<rowStart[i+1]; ++j)
                trips.push_back(Eigen::Triplet<T,sizeType>(parent._parent[i],parent._parent[value[j].second],value[j].first));
        }
        M.resize(parent._nrRow,parent._nrRow);
        M.buildFromTripletsDepulicate(trips,0.0f);
    }
    virtual AGMGLevel* getNewLevel() const {
        return new AGMGLevel();
    }
    virtual void initMem(bool redBlackTag,bool clearRB=true) {
        AMGSolver<T>::initMem(redBlackTag,clearRB);
        for(sizeType i=1; i<_nrLv; i++) {
            sizeType nrR=_levels[i]->_A.rows();
            if(!_efficient)
                ((AGMGLevel&)*(_levels[i]))._v.resize(nrR);
            ((AGMGLevel&)*(_levels[i]))._d.resize(nrR);
        }
    }
    //solve cycle
    virtual void vCycle(const FixedSparseMatrix<T,KERNEL_TYPE>& A,const Vec& r,Vec& out,const sizeType& idL) {
        if(idL < (sizeType)_nrLv-1) {
            AMGLevel& lvN=*(_levels[idL+1]);
            //pre smoothing
            for(sizeType i=0; i<(sizeType)_nrSmooth; i++)
                smooth(A,r,out,idL==0 ? _RBTag : _levels[idL]->_RBTag,true);
            //compute residue and restrict
            //compact optimized version
            computeResidueRestrict(A,r,out,lvN._r,idL+1);
            //descend
            {
                //next level v-Cycle
                KERNEL_TYPE::zero(lvN._x);
                if(_vCycle)
                    vCycle(lvN._A,lvN._r,lvN._x,idL+1);
                else if(_efficient)
                    kCycleEfficient((AGMGLevel&)lvN,lvN._x,idL+1);
                else
                    kCycle((AGMGLevel&)lvN,lvN._x,idL+1);
                //prolong add
                //compact optimized version
                computeProlongAdd(lvN._x,out,idL+1);
            }
            //post smoothing
            for(sizeType i=0; i<(sizeType)_nrSmooth; i++)
                smooth(A,r,out,idL==0 ? _RBTag : _levels[idL]->_RBTag,false);
        } else {
            //smooth for the coarsest level
            for(sizeType i=0; i<(sizeType)_nrSmoothFinal/2; i++)
                smooth(A,r,out,idL==0 ? _RBTag : _levels[idL]->_RBTag,true);
            for(sizeType i=0; i<(sizeType)_nrSmoothFinal/2; i++)
                smooth(A,r,out,idL==0 ? _RBTag : _levels[idL]->_RBTag,false);
        }
    }
    void kCycleEfficient(AGMGLevel& lv,Vec& out,const sizeType& idL) {
        if(idL == (sizeType)_nrLv-1) {
            KERNEL_TYPE::zero(out);
            //smooth for the coarsest level
            for(sizeType i=0; i<(sizeType)_nrSmoothFinal/2; i++)
                smooth(lv._A,lv._r,out,idL==0 ? _RBTag : _levels[idL]->_RBTag,true);
            for(sizeType i=0; i<(sizeType)_nrSmoothFinal/2; i++)
                smooth(lv._A,lv._r,out,idL==0 ? _RBTag : _levels[idL]->_RBTag,false);
        } else {
            for(sizeType k=0; k<_nrKCycle; k++) {
                KERNEL_TYPE::zero(out);
                vCycle(lv._A,lv._r,out,idL);
                lv._A.multiply(out,lv._d);

                T rho1,alpha1;
                if(_FCG) {
                    rho1=KERNEL_TYPE::dot(out,lv._d);
                    alpha1=KERNEL_TYPE::dot(out,lv._r);
                } else {
                    rho1=KERNEL_TYPE::dot(lv._d,lv._d);
                    alpha1=KERNEL_TYPE::dot(lv._d,lv._r);
                }

                T normR=sqrt(KERNEL_TYPE::dot(lv._r,lv._r));
                KERNEL_TYPE::addScaled(-alpha1/rho1,lv._d,lv._r);
                T normRTilde=sqrt(KERNEL_TYPE::dot(lv._r,lv._r));

                if(normRTilde < _toleranceFactor*normR) {
                    KERNEL_TYPE::scale(alpha1/rho1,out);
                } else {
                    KERNEL_TYPE::zero(lv._d);
                    vCycle(lv._A,lv._r,lv._d,idL);

                    T gamma,beta,alpha2;
                    if(_FCG) {
                        computeXtAB_XtAX(lv._A,lv._d,out,gamma,beta);
                        alpha2=KERNEL_TYPE::dot(lv._d,lv._r);
                    } else {
                        computeAXtAB_AXtAX_AXtC(lv._A,lv._d,out,lv._r,gamma,beta,alpha2);
                    }

                    T rho2=beta-gamma*gamma/rho1;
                    KERNEL_TYPE::scale(alpha1/rho1-gamma*alpha2/rho1/rho2,out);
                    KERNEL_TYPE::addScaled(alpha2/rho2,lv._d,out);
                }
            }
        }
    }
    void kCycle(AGMGLevel& lv,Vec& out,const sizeType& idL) {
        if(idL == (sizeType)_nrLv-1) {
            KERNEL_TYPE::zero(out);
            //smooth for the coarsest level
            for(sizeType i=0; i<(sizeType)_nrSmoothFinal/2; i++)
                smooth(lv._A,lv._r,out,idL==0 ? _RBTag : _levels[idL]->_RBTag,true);
            for(sizeType i=0; i<(sizeType)_nrSmoothFinal/2; i++)
                smooth(lv._A,lv._r,out,idL==0 ? _RBTag : _levels[idL]->_RBTag,false);
        } else {
            for(sizeType k=0; k<_nrKCycle; k++) {
                KERNEL_TYPE::zero(out);
                vCycle(lv._A,lv._r,out,idL);
                lv._A.multiply(out,lv._v);

                T rho1,alpha1;
                if(_FCG) {
                    rho1=KERNEL_TYPE::dot(out,lv._v);
                    alpha1=KERNEL_TYPE::dot(out,lv._r);
                } else {
                    rho1=KERNEL_TYPE::dot(lv._v,lv._v);
                    alpha1=KERNEL_TYPE::dot(lv._v,lv._r);
                }

                T normR=sqrt(KERNEL_TYPE::dot(lv._r,lv._r));
                KERNEL_TYPE::addScaled(-alpha1/rho1,lv._v,lv._r);
                T normRTilde=sqrt(KERNEL_TYPE::dot(lv._r,lv._r));

                if(normRTilde < _toleranceFactor*normR) {
                    KERNEL_TYPE::scale(alpha1/rho1,out);
                } else {
                    KERNEL_TYPE::zero(lv._d);
                    vCycle(lv._A,lv._r,lv._d,idL);

                    T gamma,beta,alpha2;
                    if(_FCG) {
                        gamma=KERNEL_TYPE::dot(lv._d,lv._v);
                        computeXtAX(lv._A,lv._d,beta);
                        alpha2=KERNEL_TYPE::dot(lv._d,lv._r);
                    } else {
                        computeAXtB_AXtAX_AXtC(lv._A,lv._d,lv._v,lv._r,gamma,beta,alpha2);
                    }

                    T rho2=beta-gamma*gamma/rho1;
                    KERNEL_TYPE::scale(alpha1/rho1-gamma*alpha2/rho1/rho2,out);
                    KERNEL_TYPE::addScaled(alpha2/rho2,lv._d,out);
                }
            }
        }
    }
    //optimized routine
    void computeXtAX(const FixedSparseMatrix<T,KERNEL_TYPE>& A,const Vec& X,T& XtAX) const {
        const std::vector<sizeType>& rowStart=A.getRowOffset();
        const FixedSparseMatrix<T,KERNEL_TYPE>::ROW& value=A.getValue();
        sizeType n=A.rows();

        T XtAXtmp=0.0f;
        sizeType i,j;
        T AX;
        #pragma omp parallel for reduction(+:XtAXtmp) private(AX,i,j)
        for(i=0; i<n; ++i) {
            AX=(T)0.0f;
            for(j=rowStart[i]; j<rowStart[i+1]; ++j)
                AX+=value[j].first*X[value[j].second];
            XtAXtmp+=AX*X[i];
        }
        XtAX=XtAXtmp;
    }
    void computeAXtB_AXtAX_AXtC(const FixedSparseMatrix<T,KERNEL_TYPE>& A,const Vec& X,const Vec& B,const Vec& C,T& AXtB,T& AXtAX,T& AXtC) const {
        const std::vector<sizeType>& rowStart=A.getRowOffset();
        const FixedSparseMatrix<T,KERNEL_TYPE>::ROW& value=A.getValue();
        sizeType n=A.rows();

        T AXtBtmp=0.0f;
        T AXtAXtmp=0.0f;
        T AXtCtmp=0.0f;
        sizeType i,j;
        T AX;
        #pragma omp parallel for reduction(+:AXtBtmp,AXtAXtmp,AXtCtmp) private(AX,i,j)
        for(i=0; i<n; ++i) {
            AX=(T)0.0f;
            for(j=rowStart[i]; j<rowStart[i+1]; ++j)
                AX+=value[j].first*X[value[j].second];
            AXtBtmp+=AX*B[i];
            AXtAXtmp+=AX*AX;
            AXtCtmp+=AX*C[i];
        }
        AXtB=AXtBtmp;
        AXtAX=AXtAXtmp;
        AXtC=AXtCtmp;
    }
    void computeXtAB_XtAX(const FixedSparseMatrix<T,KERNEL_TYPE>& A,const Vec& X,const Vec& B,T& XtAB,T& XtAX) const {
        const std::vector<sizeType>& rowStart=A.getRowOffset();
        const FixedSparseMatrix<T,KERNEL_TYPE>::ROW& value=A.getValue();
        sizeType n=A.rows();

        T XtABtmp=0.0f;
        T XtAXtmp=0.0f;
        sizeType i,j,col;
        T AX,AB,val;
        #pragma omp parallel for reduction(+:XtABtmp,XtAXtmp) private(AX,AB,i,j,col,val)
        for(i=0; i<n; ++i) {
            AX=(T)0.0f;
            AB=(T)0.0f;
            for(j=rowStart[i]; j<rowStart[i+1]; ++j) {
                val=value[j].first;
                col=value[j].second;
                AX+=val*X[col];
                AB+=val*B[col];
            }
            XtABtmp+=AB*X[i];
            XtAXtmp+=AX*X[i];
        }
        XtAB=XtABtmp;
        XtAX=XtAXtmp;
    }
    void computeAXtAB_AXtAX_AXtC(const FixedSparseMatrix<T,KERNEL_TYPE>& A,const Vec& X,const Vec& B,const Vec& C,T& AXtAB,T& AXtAX,T& AXtC) const {
        const std::vector<sizeType>& rowStart=A.getRowOffset();
        const FixedSparseMatrix<T,KERNEL_TYPE>::ROW& value=A.getValue();
        sizeType n=A.rows();

        T AXtABtmp=0.0f;
        T AXtAXtmp=0.0f;
        T AXtCtmp=0.0f;
        T AX,AB,val;
        sizeType i,j,col;
        #pragma omp parallel for reduction(+:AXtABtmp,AXtAXtmp,AXtCtmp) private(AX,AB,i,j,col,val)
        for(i=0; i<n; ++i) {
            AX=(T)0.0f;
            AB=(T)0.0f;
            for(j=rowStart[i]; j<rowStart[i+1]; ++j) {
                val=value[j].first;
                col=value[j].second;
                AX+=val*X[col];
                AB+=val*B[col];
            }
            AXtABtmp+=AX*AB;
            AXtAXtmp+=AX*AX;
            AXtCtmp+=AX*C[i];
        }
        AXtAB=AXtABtmp;
        AXtAX=AXtAXtmp;
        AXtC=AXtCtmp;
    }
    void computeResidueRestrict(const FixedSparseMatrix<T,KERNEL_TYPE>& A,const Vec& R,const Vec& X,Vec& XN,const sizeType& lv) const {
        const std::vector<sizeType>& rowStart=A.getRowOffset();
        const FixedSparseMatrix<T,KERNEL_TYPE>::ROW& value=A.getValue();

        const std::vector<sizeType>& parent=((AGMGLevel&)*(_levels[lv]))._parent._parent;
        const sizeType nr=(sizeType)parent.size();

        KERNEL_TYPE::zero(XN);
        sizeType i,j,p;
        T tmp;
        #pragma omp parallel for private(i,j,p,tmp)
        for(i=0; i<nr; i++) {
            p=parent[i];
            if(p >= 0) {
                tmp=R[i];
                for(j=rowStart[i]; j<rowStart[i+1]; ++j)
                    tmp-=value[j].first*X[value[j].second];
                OMP_ATOMIC_
                XN[p]+=tmp;
            }
        }
    }
    void computeRestrict(const Vec& X,Vec& XN,const sizeType& lv) const {
        const std::vector<sizeType>& parent=((AGMGLevel&)*(_levels[lv]))._parent._parent;
        const sizeType nr=(sizeType)parent.size();

        KERNEL_TYPE::zero(XN);
        sizeType i,p;
        #pragma omp parallel for private(i,p)
        for(i=0; i<nr; i++) {
            p=parent[i];
            if(p >= 0) {
                OMP_ATOMIC_
                XN[p]+=X[i];
            }
        }
    }
    void computeProlongAdd(const Vec& XN,Vec& X,const sizeType& lv) const {
        const std::vector<sizeType>& parent=((AGMGLevel&)*(_levels[lv]))._parent._parent;
        const sizeType nr=(sizeType)parent.size();
        sizeType i,p;
        #pragma omp parallel for private(i,p)
        for(i=0; i<nr; i++) {
            p=parent[i];
            if(p >= 0)
                X[i]+=XN[p];
        }
    }
public:
    //for general debug
    template <typename T2>
    void debugLevel(const Grid<sizeType,T2>& off) const {
        ASSERT(off.getDim() == 2);
        std::map<sizeType,Vec3i> offMap;
        const Vec3i nrPoint=off.getNrPoint();
        for(sizeType x=0; x<nrPoint.x(); x++)
            for(sizeType y=0; y<nrPoint.y(); y++) {
                Vec3i id(x,y,0);
                offMap[off.get(id)]=id;
            }

        std::vector<sizeType> parent;
        parent.resize(_A->rows());
        for(sizeType i=0; i<(sizeType)parent.size(); i++)
            parent[i]=i;

        for(sizeType i=0; i<_nrLv; i++) {
            //write VTK
            debugWriteLevel(offMap,off,parent,i,(i == 0) ? _RBTag : _levels[i]->_RBTag);
            //to next level
            if(i < (sizeType)_nrLv-1) {
                const std::vector<sizeType>& pp=((AGMGLevel&)(*(_levels[i+1])))._parent._parent;
                for(sizeType j=0; j<(sizeType)parent.size(); j++)
                    if(parent[j]>=0)
                        parent[j]=pp[parent[j]];
            }
        }
    }
    template <typename T2>
    void debugWriteLevel(const std::map<sizeType,Vec3i>& offMap,const Grid<sizeType,T2>& off,const std::vector<sizeType>& parent,const sizeType& lv,const std::vector<bool>& RBTag) const {
        Grid<sizeType,T2> tag=off;
        tag.init(-1);
        for(sizeType i=0; i<(sizeType)parent.size(); i++) {
            std::map<sizeType,Vec3i>::const_iterator iter=offMap.find(i);
            if(iter != offMap.end())
                tag.get(iter->second)=parent[i];
        }

        ostringstream oss;
        oss << "./level" << lv << ".vtk";
        VTKWriter<T2> os("AGMG Level",oss.str(),true);

        std::vector<ScalarUtil<T2>::ScalarVec3> voxels;
        std::vector<T2> color;
        const Vec3i nrPoint=off.getNrPoint();
        for(sizeType x=0; x<nrPoint.x(); x++)
            for(sizeType y=0; y<nrPoint.y(); y++) {
                Vec3i id(x,y,0);
                if(off.get(id) == -1 || tag.get(id) == -1)
                    color.push_back(0.0f);
                else if(!RBTag.empty() && RBTag[tag.get(id)])
                    color.push_back(0.5f);
                else color.push_back(1.0f);

                ScalarUtil<T2>::ScalarVec3 minOff=off.getCellSize()*0.5f;
                ScalarUtil<T2>::ScalarVec3 maxOff=off.getCellSize()*0.5f;
                if(tag.isSafeIndex(Vec3i(x-1,y,0)) && tag.get(Vec3i(x-1,y,0)) != tag.get(id))
                    minOff.x()*=0.8f;
                if(tag.isSafeIndex(Vec3i(x+1,y,0)) && tag.get(Vec3i(x+1,y,0)) != tag.get(id))
                    maxOff.x()*=0.8f;
                if(tag.isSafeIndex(Vec3i(x,y-1,0)) && tag.get(Vec3i(x,y-1,0)) != tag.get(id))
                    minOff.y()*=0.8f;
                if(tag.isSafeIndex(Vec3i(x,y+1,0)) && tag.get(Vec3i(x,y+1,0)) != tag.get(id))
                    maxOff.y()*=0.8f;

                ScalarUtil<T2>::ScalarVec3 minC=off.getPt(id)-minOff;
                ScalarUtil<T2>::ScalarVec3 maxC=off.getPt(id)+maxOff;
                minC.z()=-off.getCellSize().maxCoeff()*0.01f;
                maxC.z()=+off.getCellSize().maxCoeff()*0.01f;

                voxels.push_back(minC);
                voxels.push_back(maxC);
            }
        os.appendVoxels(voxels.begin(),voxels.end(),false);
        os.appendCustomData("color",color.begin(),color.end());
    }
    //for profile
    void printMemoryConsumption() const {
        INFOV("Memory Consumption Multigrid")
        T totalMBA=0.0f,totalMBXRVD=0.0f;
        for(sizeType i=0; i<(sizeType)_nrLv; i++) {
            T MBA,MBXRVD;
            ((AGMGLevel&)*(_levels[i])).getMemoryConsumption(MBA,MBXRVD);
            INFOV("Level %d: %fMB(%fGB) for A and %fMB(%fGB) for X,R,V,B",i,MBA,MBA/1024.0f,MBXRVD,MBXRVD/1024.0f)
            totalMBA+=MBA;
            totalMBXRVD+=MBXRVD;
        }
        INFOV("Total: %fMB(%fGB) for A and %fMB(%fGB) for X,R,V,B",totalMBA,totalMBA/1024.0f,totalMBXRVD,totalMBXRVD/1024.0f)
    }
protected:
    //data
    sizeType _nrKCycle;
    bool _FCG;
    bool _efficient;
    bool _checkDD;
    bool _vCycle;
	bool _reuseAggr;
};

PRJ_END

#endif