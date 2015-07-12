#ifndef AMG_SOLVER_H
#define AMG_SOLVER_H

#include "LinearSolver.h"
#include "DisjointSet.h"
#include "Heap.h"
#include "GridOp.h"
#include <queue>

PRJ_BEGIN

//An classical AMG solver
template <typename T>
class AMGSolver : public Solver<T,Kernel<T> >
{
public:
    typedef Kernel<T> KERNEL_TYPE;
    using typename Solver<T,KERNEL_TYPE>::Vec;
    using typename Solver<T,KERNEL_TYPE>::SOLVER_RESULT;
    using Solver<T,KERNEL_TYPE>::_cb;
    using Solver<T,KERNEL_TYPE>::_iterationsOut;
    using Solver<T,KERNEL_TYPE>::SUCCESSFUL;
    using Solver<T,KERNEL_TYPE>::NOT_CONVERGENT;
	struct AMGLevel
	{
		FixedSparseMatrix<T,KERNEL_TYPE> _IT;
		FixedSparseMatrix<T,KERNEL_TYPE> _A;
		std::vector<sizeType> _idC;
        std::vector<bool> _RBTag;
		Vec _x;
		Vec _r;
	};
	AMGSolver()
	{
        _nrLv=0;
		_alpha=0.2f;
		_nrSmooth=2;
		_nrSmoothFinal=10;
        _asPrecon=true;

		setSolverParameters(1e-5f,1000);
	}
    virtual ~AMGSolver(){}
	virtual void setSolverParameters(T toleranceFactor,sizeType maxIterations)
	{
		_toleranceFactor=toleranceFactor;
		if(_toleranceFactor<1e-30f) 
			_toleranceFactor=1e-30f;
		_maxIterations=maxIterations;
	}
	virtual void setMatrix(const SparseMatrix<T,KERNEL_TYPE>& matrix,bool syncPrecon=true)
	{
        _AStorage.constructFromMatrix(matrix);
		setMatrix(_AStorage,syncPrecon);
	}
	virtual void setMatrix(const FixedSparseMatrix<T,KERNEL_TYPE>& matrix,bool syncPrecon=true)
	{
		FixedSparseMatrix<T,KERNEL_TYPE> currentMatrix=matrix;
		boost::shared_ptr<AMGLevel> currentLevel;
		std::vector<T> maxCouple;
		std::vector<sizeType> idC;

		//build base level
        _nrLv=0;
		currentLevel.reset((AMGLevel*)NULL);
        _A=&matrix;

        if((sizeType)_levels.size() <= _nrLv)
            _levels.resize(_nrLv+1);
		_levels[_nrLv]=currentLevel;
		_nrLv++;

		//coarsen
		sizeType lastNrC=currentMatrix.rows();
		while(true)
		{
			sizeType nrC;
			TBEG("Coarsen")
			nrC=coarsen(currentMatrix,maxCouple,idC);
			TEND
			if(nrC <= 1 || nrC >= lastNrC*0.75f)
				break;
			else{
				TBEG("Build Level")
				currentLevel.reset(getNewLevel());
				buildMatrix(*currentLevel,currentMatrix,maxCouple,idC,nrC);
				//currentLevel->init();
				TEND

                if((sizeType)_levels.size() <= _nrLv)
                    _levels.resize(_nrLv+1);
		        _levels[_nrLv]=currentLevel;
				_nrLv++;
			}
			lastNrC=nrC;
		}
        initMem(true);
	}
	virtual void setMatrix(const std::vector<FixedSparseMatrix<T,KERNEL_TYPE> >& As,const std::vector<FixedSparseMatrix<T,KERNEL_TYPE> >& ITs)
	{
        _nrLv=As.size();
        if((sizeType)_levels.size() < _nrLv)
		    _levels.resize(As.size());
        _levels[0].reset((AMGLevel*)NULL);
        _A=&(As[0]);

		for(sizeType i=1;i<(sizeType)As.size();i++)
		{
			_levels[i].reset(getNewLevel());
			_levels[i]->_A=As[i];
			if(i > 0)
				_levels[i]->_IT=ITs[i-1];
		}
        initMem(true);
	}
    //you just provide a aggregation matrix
	//the final prolongation matrix=(ID-(2.0f/3.0f)*D^-1*A)*I
	//the final restriction matrix=IT*(ID-(2.0f/3.0f)*A^T*D^-1)
	virtual void setMatrix(const FixedSparseMatrix<T,KERNEL_TYPE>& A0,std::vector<FixedSparseMatrix<T,KERNEL_TYPE> >& ITs,bool smooth=true)
	{
		FixedSparseMatrix<T,KERNEL_TYPE> lastA=A0;
        _nrLv=(sizeType)ITs.size()+1;
		if((sizeType)_levels.size() < _nrLv)
            _levels.resize(_nrLv);
		//base level
		_levels[0].reset((AMGLevel*)NULL);
        _A=&A0;
		
		//build coarser levels
		for(sizeType i=0;i<(sizeType)ITs.size();i++)
		{
			_levels[i+1].reset(getNewLevel());

			//build the coarse level restriction matrix
            if(smooth)
			{
				FixedSparseMatrix<T,KERNEL_TYPE> ITA,WInvD,ITAInvD;
				ITs[i].multiplyMatrixTransposed(lastA,ITA);
				{
					WInvD.resize(lastA.rows());
					for(sizeType r=0;r<lastA.rows();r++)
					{
						T D=lastA(r,r);
						if(D != 0.0f)
							D=(2.0f/3.0f)/D;
						WInvD.setElement(r,r,D);
					}
				}
				ITA.multiplyMatrix(WInvD,ITAInvD);
				ITs[i].sub(ITAInvD);
			}

			//build the coarse level matrix
			FixedSparseMatrix<T,KERNEL_TYPE> ITA;
			ITs[i].multiplyMatrix(lastA,ITA);
			ITA.multiplyMatrixTransposed(ITs[i],lastA);
			_levels[i+1]->_A=lastA;
			//build the rest
			_levels[i+1]->_IT=ITs[i];
		}
        initMem(true);
	}
	virtual SOLVER_RESULT solve(const Vec& rhs,Vec& result)
	{
        if(_asPrecon)
        {
            KERNEL_TYPE::zero(result);
            vCycle(*_A,rhs,result,0);
            return SUCCESSFUL;
        }

		if(_cb)
			_cb->reset();

		result.resize(_A->rows());
		KERNEL_TYPE::zero(result);
		sizeType iteration;
		for(iteration=0; iteration<_maxIterations; ++iteration)
		{
			//check convergency
			_residualOut=findMaximalResidue(*_A,rhs,result);
			if(_cb)
				(*_cb)(result,_residualOut,iteration);
			if(_residualOut < _toleranceFactor)
			{
				_iterationsOut=iteration;
				return SUCCESSFUL;
			}

			//v-Cycle
			vCycle(*_A,rhs,result,0);
		}

		_iterationsOut=iteration;
		return NOT_CONVERGENT;
	}
	sizeType n() const{return _levels[0]->_A.row();}
    void setAsPreconditioner(bool asPrecon){_asPrecon=asPrecon;}
    sizeType nrLv() const{return _nrLv;}
	void setNrSmooth(sizeType nr){_nrSmooth=nr;}
	void setNrSmoothFinal(sizeType nr){_nrSmoothFinal=nr;}
	void debugFillin() const
	{
		if(_A)_A->printFillin();
		for(sizeType i=0;i<(sizeType)_levels.size();i++)
			if(_levels[i])_levels[i]->_A.printFillin();
	}
protected:
	virtual sizeType coarsen(const FixedSparseMatrix<T,KERNEL_TYPE>& sm,std::vector<T>& maxCouple,std::vector<sizeType>& idC) const
	{
#define FAR 0
#define COARSE 1
#define FINE 2
		std::vector<sizeType> heap;
		std::vector<sizeType> heapOff(sm.rows());
		std::vector<sizeType> lambda(sm.rows(),0);
		std::vector<std::set<sizeType> > sT(sm.rows());
		maxCouple.resize(sm.rows());
		idC.assign(sm.rows(),FAR);

		//coarsen phase I
		for(sizeType i=0;i<sm.rows();i++)
		{
			T maxVal=0.0f;
			for(ConstSMIterator<T> beg=sm.begin(i),end=sm.end(i);beg!=end;++beg)
				if(beg.col() != i)
					maxVal=std::max(maxVal,-(*beg));
			maxCouple[i]=maxVal;

			for(ConstSMIterator<T> beg=sm.begin(i),end=sm.end(i);beg!=end;++beg)
				if(beg.col() != i && -(*beg) >= maxCouple[i]*_alpha)
				{
					lambda[beg.col()]++;
					sT[beg.col()].insert(i);
				}
		}
		for(sizeType i=0;i<sm.rows();i++)
			pushHeapMaxDef(lambda,heapOff,heap,i);
		while(!heap.empty())
		{
			sizeType err;
			sizeType i=popHeapMaxDef(lambda,heapOff,heap,err);
			if(idC[i] == FAR)
			{
				if(maxCouple[i] == 0.0f)
				{
					idC[i]=FINE;
					continue;
				}

				idC[i]=COARSE;
				for(std::set<sizeType>::const_iterator beg=sT[i].begin(),end=sT[i].end();beg!=end;beg++)
				{
					sizeType j=*beg;
					if(idC[j] == FAR)
					{
						idC[j]=FINE;
						for(ConstSMIterator<T> beg=sm.begin(j),end=sm.end(j);beg!=end;++beg)
							if(beg.col() != j && -(*beg) >= maxCouple[j]*_alpha)
							{
								lambda[beg.col()]++;
								updateHeapMaxDef(lambda,heapOff,heap,beg.col());
							}
					}
				}
			}
		}

		//coarsen phase II
		for(sizeType i=0;i<sm.rows();i++)
		{
			if(idC[i] == FINE)
			{
				if(maxCouple[i] == 0.0f)
					continue;

				std::set<sizeType> coarse;
				std::vector<sizeType> fine;

				//find fine and coarse dependent set
				{
					for(ConstSMIterator<T> beg=sm.begin(i),end=sm.end(i);beg!=end;++beg)
						if(beg.col() != i && -(*beg) >= maxCouple[i]*_alpha)
						{
							if(idC[beg.col()] == COARSE)
								coarse.insert(beg.col());
							else fine.push_back(beg.col());
						}
				}

				for(sizeType j=0;j<(sizeType)fine.size();j++)
				{
					sizeType idF=fine[j];
					if(maxCouple[idF] == 0.0f)
						continue;

					bool common=false;
					for(ConstSMIterator<T> beg=sm.begin(idF),end=sm.end(idF);beg!=end;++beg)
						if(beg.col() != idF && -(*beg) >= maxCouple[idF]*_alpha)
						{
							if(coarse.find(beg.col()) != coarse.end())
							{
								common=true;
								break;
							}
						}
					if(!common)
					{
						idC[idF]=COARSE;
						coarse.insert(idF);
						fine.erase(fine.begin()+j);
						j--;
					}
				}
			}
		}

		//build the coarse index in interpolation and restriction matrix
		sizeType j=0;
		for(sizeType i=0;i<sm.rows();i++)
		{
			if(idC[i] == COARSE)
				idC[i]=j++;
			else idC[i]=-1;
		}
		return j;
	}
	virtual void buildMatrix(AMGLevel& level,FixedSparseMatrix<T,KERNEL_TYPE>& A,const std::vector<T>& maxCouple,const std::vector<sizeType>& idC,const sizeType& nrC) const
	{
		level._idC=idC;

		//build restriction matrix
		FixedSparseMatrix<T,KERNEL_TYPE> IT;
		IT.resize(nrC,A.rows());
		vector<Eigen::Triplet<scalarD,sizeType> > trips;
		trips.reserve(A.getValue().size());
		for(sizeType i=0;i<A.rows();i++)
		{
			if(idC[i] >= 0)
			{
				trips.push_back(Eigen::Triplet<scalarD,sizeType>(idC[i],i,1.0f));
			}
			else
			{
				if(maxCouple[i] == 0.0f)
					continue;

				for(SMIterator<T> beg=A.begin(i),end=A.end(i);beg!=end;++beg)
				{
					sizeType idJ=beg.col();
					sizeType idCJ=idC[idJ];
					if(idCJ >= 0)
					{
						T w=weight(i,idJ,A,maxCouple,idC);
						trips.push_back(Eigen::Triplet<scalarD,sizeType>(idCJ,i,w));
					}
				}
			}
		}IT.buildFromTripletsDepulicate(trips,0.0f);

		//set A=I^T*A*I
		FixedSparseMatrix<T,KERNEL_TYPE> ITA;
		IT.multiplyMatrix(A,ITA);
		ITA.multiplyMatrixTransposed(IT,A);

		//set to level
		level._A=A;
		level._IT=IT;
	}
	virtual T weight(const sizeType& i,const sizeType& j,const FixedSparseMatrix<T,KERNEL_TYPE>& A,const std::vector<T>& maxCouple,const std::vector<sizeType>& idC) const
	{
		std::vector<sizeType> strong;
		std::vector<sizeType> coarse;
		std::vector<sizeType> weak;
		for(ConstSMIterator<T> beg=A.begin(i),end=A.end(i);beg!=end;++beg)
		{
			if(beg.col() != i)
			{
				if(idC[beg.col()] >= 0)
					coarse.push_back(beg.col());
				else if(-(*beg) >= maxCouple[i]*_alpha)
					strong.push_back(beg.col());
				else weak.push_back(beg.col());
			}
		}

		T wWeak=A(i,i);
		for(sizeType k=0;k<(sizeType)weak.size();k++)
			wWeak+=A(i,weak[k]);

		T wStrong=0.0f;
		for(sizeType k=0;k<(sizeType)strong.size();k++)
		{
			T wStrongCoarse=0.0f;
			for(sizeType l=0;l<(sizeType)coarse.size();l++)
				wStrongCoarse+=A(strong[k],coarse[l]);

			//this assert should always hold because we work on an M matrix
			//and our coarsening algorithm ensure that For each point i \in F,
			//the strong dependent set S_i and each fine point j \in S_i:
			//if j \in F, than j must strongly depnds on (C \intersect S_i)
			ASSERT(std::abs(wStrongCoarse) > 0.0f)
			wStrong+=A(i,strong[k])*A(strong[k],j)/wStrongCoarse;
		}
		return -(A(i,j)+wStrong)/wWeak;
	}
    virtual void vCycle(const FixedSparseMatrix<T,KERNEL_TYPE>& A,const Vec& r,Vec& out,const sizeType& idL)
	{
        if(_xBK.size() < A.rows())
            _xBK.resize(A.rows());

		if(idL < (sizeType)_nrLv-1)
		{
			AMGLevel& lvN=*(_levels[idL+1]);
			//pre smoothing
			for(sizeType i=0;i<(sizeType)_nrSmooth;i++)
				smooth(A,r,out,idL==0 ? _RBTag : _levels[idL]->_RBTag,true);
			//compute residue
            Eigen::Block<Vec> xBK=_xBK.block(0,0,A.rows(),1);
			xBK=r;
			A.multiplySubtract(out,xBK);
			//descend
			{
				//restrict
				lvN._IT.multiply(xBK,lvN._r);
				//next level v-Cycle
				KERNEL_TYPE::zero(lvN._x);
				vCycle(lvN._A,lvN._r,lvN._x,idL+1);
				//prolong add
				lvN._IT.multiplyTransposeAdd(lvN._x,out);
			}
			//post smoothing
			for(sizeType i=0;i<(sizeType)_nrSmooth;i++)
				smooth(A,r,out,idL==0 ? _RBTag : _levels[idL]->_RBTag,false);
		}
		else
		{
			//smooth for the coarsest level
			for(sizeType i=0;i<(sizeType)_nrSmoothFinal/2;i++)
				smooth(A,r,out,idL==0 ? _RBTag : _levels[idL]->_RBTag,true);
            for(sizeType i=0;i<(sizeType)_nrSmoothFinal/2;i++)
				smooth(A,r,out,idL==0 ? _RBTag : _levels[idL]->_RBTag,false);
		}
	}
    virtual AMGLevel* getNewLevel() const{return new AMGLevel();}
    virtual void initMem(bool redBlackTag,bool clearRB=true)
    {
        if(clearRB)
            _RBTag.clear();
        if(redBlackTag)
            decideRB(*_A,_RBTag);
        for(sizeType i=1;i<_nrLv;i++)
        {
            sizeType nrR=_levels[i]->_A.rows();
            _levels[i]->_x.resize(nrR);
            _levels[i]->_r.resize(nrR);
            if(clearRB)
                _levels[i]->_RBTag.clear();
            if(redBlackTag)
                decideRB(_levels[i]->_A,_levels[i]->_RBTag);
        }
    }
    virtual void decideRB(const FixedSparseMatrix<T,Kernel<T> >& A,std::vector<bool>& tag) const
    {
        const std::vector<sizeType>& rowStart=A.getRowOffset();
        const FixedSparseMatrix<T,KERNEL_TYPE>::ROW& value=A.getValue();
        sizeType n=A.rows();

        sizeType i,j,k,C,sz,node;
        unsigned char curr;
        std::queue<sizeType> q;
        std::vector<unsigned char> tagTmp(n,0);
        tag.resize(n);
        for(i=0;i<n;i++)
        {
            //an untagged tree in the forest
            if(tagTmp[i] == 0)
            {
                //tag this tree by BFS
                tagTmp[i]=1;
                curr=2;
                q.push(i);
                while(!q.empty())
                {
                    sz=(sizeType)q.size();
                    for(j=0;j<sz;j++)
                    {
                        node=q.front();q.pop();
                        //tag everyone connected
                        for(k=rowStart[node];k<rowStart[node+1];k++)
                        {
                            C=value[k].second;
                            if(tagTmp[C] == 0){
                                tagTmp[C]=curr;
                                q.push(C);
                            }
                        }
                    }
                    curr=3-curr;
                }
            }
            tag[i]=tagTmp[i]==1;
        }
    }
    //optimized routines
	virtual void smooth(const FixedSparseMatrix<T,KERNEL_TYPE>& A,const Vec& R,Vec& X,const std::vector<bool>& RBTag,bool positive) const
	{
        if(!RBTag.empty()){
            RBGS(A,R,X,RBTag,positive,A.rows());
            RBGS(A,R,X,RBTag,!positive,A.rows());
        }else{
            gaussSeidel(A,R,X,positive,A.rows());
		    dampedJacobi(A,R,X,2.0f/3.0f,A.rows());
        }
	}
    void RBGS(const FixedSparseMatrix<T,KERNEL_TYPE>& A,const Vec& R,Vec& X,const std::vector<bool>& RBTag,bool red,sizeType nrRow) const
	{
        const std::vector<sizeType>& rowStart=A.getRowOffset();
        const FixedSparseMatrix<T,KERNEL_TYPE>::ROW& value=A.getValue();

        sizeType i,j,C;
        T D;
        #pragma omp parallel for private(i,j,C,D)
		for(i=0;i<nrRow;i++)
		{
            if(RBTag[i] != red)
                continue;

			D=0.0f;
			X[i]=R[i];
			for(j=rowStart[i]; j<rowStart[i+1]; ++j){
                C=value[j].second;
				if(C == i) D=value[j].first;
				else X[i]-=value[j].first*X[C];
			}
			if(D > 0.0f)X[i]/=D;
			else X[i]=0.0f;
		}
	}
    void gaussSeidel(const FixedSparseMatrix<T,KERNEL_TYPE>& A,const Vec& R,Vec& X,bool positive,sizeType nrRow) const
	{
        const std::vector<sizeType>& rowStart=A.getRowOffset();
        const FixedSparseMatrix<T,KERNEL_TYPE>::ROW& value=A.getValue();

        sizeType from,to,delta;
        if(!positive){
            from=0;
            to=nrRow-1;
            delta=1;
        }else{
            from=nrRow-1;
            to=0;
            delta=-1;
        }

        sizeType i,j,C;
        T D;
		for(i=from;i!=to;i+=delta)
		{
			D=0.0f;
			X[i]=R[i];
			for(j=rowStart[i]; j<rowStart[i+1]; ++j){
                C=value[j].second;
				if(C == i) D=value[j].first;
				else X[i]-=value[j].first*X[C];
			}
			if(D > 0.0f)X[i]/=D;
			else X[i]=0.0f;
		}
	}
	void dampedJacobi(const FixedSparseMatrix<T,KERNEL_TYPE>& A,const Vec& R,Vec& X,const T& w,sizeType nrRow) const
	{
        Vec xBK;
        if(xBK.size() != X.size())
            xBK.resize(X.size());

        const std::vector<sizeType>& rowStart=A.getRowOffset();
        const FixedSparseMatrix<T,KERNEL_TYPE>::ROW& value=A.getValue();

        sizeType i,j,C;
        T D;
        //#pragma omp parallel for private(i,j,C,D)
		for(i=0;i<nrRow;i++)
		{
			D=0.0f;
			xBK[i]=R[i];
			for(j=rowStart[i]; j<rowStart[i+1]; ++j){
                C=value[j].second;
				if(C == i) 
                    D=value[j].first;
				xBK[i]-=value[j].first*X[C];
			}
			if(D > 0.0f) xBK[i]*=(w/D);
			else xBK[i]=0.0f;
		}
		KERNEL_TYPE::add(xBK,X,X);
	}
    T findMaximalResidue(const FixedSparseMatrix<T,KERNEL_TYPE>& A,const Vec& R,const Vec& X) const
    {
        const std::vector<sizeType>& rowStart=A.getRowOffset();
        const FixedSparseMatrix<T,KERNEL_TYPE>::ROW& value=A.getValue();
        sizeType n=A.rows();

        T maxResShared=0.0f;
        T maxRes=0.0f;
        #pragma omp parallel shared(maxResShared) firstprivate(maxRes)
        {
            #pragma omp for nowait
            for(sizeType i=0;i<n;i++) {
                T res=(T)R[i];
                for(sizeType j=rowStart[i]; j<rowStart[i+1]; ++j)
                    res-=value[j].first*X[value[j].second];
                maxRes=std::max(std::abs(res),maxRes);
            }
            OMP_CRITICAL_
            {
                maxResShared=std::max(maxRes,maxResShared);
            }
        }
        return maxResShared;
    }
public:
	//for grid problem debug
	void debug(const Vec3i& nrCell) const
	{
        typedef typename ScalarUtil<T>::ScalarVec3 Vec3Type;
		Grid<T,T> cells;
        cells.reset(nrCell,BBox<T>(Vec3Type(0.0f,0.0f,0.0f),Vec3Type((T)nrCell.x(),(T)nrCell.y(),0.0f)),1.0f);

		for(sizeType i=1;i<(sizeType)_nrLv;i++)
		{
			debugLevel(cells,*(_levels[i]),nrCell);
			ostringstream oss;oss << "./level" << i << ".vtk";
			GridOp<T,T>::write2DScalarBarChartVTK(oss.str(),cells,false);
		}
		debugBaseLevel(cells,*(_levels[0]),nrCell);
		GridOp<T,T>::write2DScalarBarChartVTK("./level0.vtk",cells,false);
	}
	void debugBaseLevel(Grid<T,T>& cells,const AMGLevel& base,const Vec3i& nrCell) const
	{
		for(sizeType x=0;x<nrCell.x();x++)
		for(sizeType y=0;y<nrCell.y();y++)
		{
			const sizeType id=x*nrCell.y()+y;
			if(base._A.numElement(id) > 0)
				cells.get(Vec3i(x,y,0))=1.0f;
			else cells.get(Vec3i(x,y,0))=0.0f;
		}
	}
	void debugLevel(Grid<T,T>& cells,const AMGLevel& base,const Vec3i& nrCell) const
	{
		sizeType idC=0;
		for(sizeType x=0;x<nrCell.x();x++)
		for(sizeType y=0;y<nrCell.y();y++)
		{
			if(cells.get(Vec3i(x,y,0)) == 1.0f)
			{
				if(base._idC[idC] < 0)
					cells.get(Vec3i(x,y,0))=0.0f;
				idC++;
			}
		}
	}
    static void debugRBTag(const FixedSparseMatrix<T,Kernel<T> >& A,const std::vector<bool>& RBTag)
    {
        for(sizeType i=0;i<A.rows();i++)
        {
            bool currRB=RBTag[i];
            for(ConstSMIterator<T> beg=A.begin(i),end=A.end(i);beg!=end;++beg)
                if(beg.col() != i)
                    ASSERT(currRB != RBTag[beg.col()]);
        }
    }
    //for debug case:
	//mask == 0 means neumann boundary
	//mask == 1 means variable
	//mask == 2 meaks dirichlet boundary
	void buildTemplateMatrix(const Vec3i& size,const Grid<T,T>& tpl,const Vec3i& pos,const Grid<T,T>* mask)
	{
		FixedSparseMatrix<T,KERNEL_TYPE> mat;
		mat.resize(size.x()*size.y());
		for(sizeType x=0;x<size.x();x++)
		for(sizeType y=0;y<size.y();y++)
		{
#define GI(X,Y) ((X)*size.y()+(Y))
			if(!mask || mask->get(Vec3i(x,y,0)) == 1.0f)
			{
				T ctrVal=tpl.get(Vec3i(pos.x(),pos.y(),0));
				for(sizeType xx=0;xx<tpl.getNrPoint().x();xx++)
				for(sizeType yy=0;yy<tpl.getNrPoint().y();yy++)
				{
					Vec3i pt(x+xx-pos.x(),y+yy-pos.y(),0);
					if(pt.x() >= 0 && pt.x() < size.x() && pt.y() >= 0 && pt.y() < size.y())
					{
						if(!mask || mask->get(pt) == 1.0f)
							mat.setElement(GI(x,y),GI(pt.x(),pt.y()),tpl.get(Vec3i(xx,yy,0)));
						else if(mask->get(pt) == 0.0f)
							ctrVal-=std::abs(tpl.get(Vec3i(xx,yy,0)));
					}else ctrVal-=std::abs(tpl.get(Vec3i(xx,yy,0)));
				}
				mat.setElement(GI(x,y),GI(x,y),ctrVal);
			}
		}
		setMatrix(mat);
	}
	void buildUniformLaplaceMatrix4(const Vec3i& size,const Grid<T,T>* mask)
	{
		Grid<T,T> tpl;
		typedef typename ScalarUtil<T>::ScalarVec3 Vec3Type;
        tpl.reset(Vec3i(3,3,0),BBox<T>(Vec3Type(0.0f,0.0f,0.0f),Vec3Type(1.0f,1.0f,0.0f)),0.0f);
		tpl.get(Vec3i(0,0,0))= 0.0f;
		tpl.get(Vec3i(1,0,0))=-1.0f;
		tpl.get(Vec3i(2,0,0))= 0.0f;
		tpl.get(Vec3i(0,1,0))=-1.0f;
		tpl.get(Vec3i(1,1,0))= 4.0f;
		tpl.get(Vec3i(2,1,0))=-1.0f;
		tpl.get(Vec3i(0,2,0))= 0.0f;
		tpl.get(Vec3i(1,2,0))=-1.0f;
		tpl.get(Vec3i(2,2,0))= 0.0f;
		buildTemplateMatrix(size,tpl,Vec3i(1,1,0),mask);
	}
	void buildUniformLaplaceMatrix8(const Vec3i& size,const Grid<T,T>* mask)
	{
        Grid<T,T> tpl;
		typedef typename ScalarUtil<T>::ScalarVec3 Vec3Type;
        tpl.reset(Vec3i(3,3,0),BBox<T>(Vec3Type(0.0f,0.0f,0.0f),Vec3Type(1.0f,1.0f,0.0f)),0.0f);
		tpl.get(Vec3i(0,0,0))=-0.25f;
		tpl.get(Vec3i(1,0,0))=-0.5f;
		tpl.get(Vec3i(2,0,0))=-0.25f;
		tpl.get(Vec3i(0,1,0))=-0.5f;
		tpl.get(Vec3i(1,1,0))= 3.0f;
		tpl.get(Vec3i(2,1,0))=-0.5f;
		tpl.get(Vec3i(0,2,0))=-0.25f;
		tpl.get(Vec3i(1,2,0))=-0.5f;
		tpl.get(Vec3i(2,2,0))=-0.25f;
		buildTemplateMatrix(size,tpl,Vec3i(1,1,0),mask);
	}
    void profileVCycle()
    {
        BEGIN_PROFILE
        for(sizeType i=0;i<10;i++)
            vCycle(*(_levels[0]),_levels[0]->_x,0);
        END_PROFILE("100 VCycle Time")
    }
protected:
	//AMG structure
    const FixedSparseMatrix<T,Kernel<T> > *_A;
    FixedSparseMatrix<T,Kernel<T> > _AStorage;
    std::vector<bool> _RBTag;
    //deeper level
	std::vector<boost::shared_ptr<AMGLevel> > _levels;
    sizeType _nrLv;
    //multigrid options
	T _alpha;
	sizeType _nrSmooth;
	sizeType _nrSmoothFinal;
    bool _asPrecon;
    //for termination
	T _residualOut;
	T _toleranceFactor;
	sizeType _maxIterations;
    //helper
    Vec _xBK;
};

PRJ_END

#endif
