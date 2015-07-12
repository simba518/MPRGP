#ifndef AINV_PRECONDITIONER_H
#define AINV_PRECONDITIONER_H

#include "LinearSolver.h"

PRJ_BEGIN

template <typename T,typename KERNEL_TYPE=Kernel<T> >
struct SymAINVPreconSolver : public Solver<T,KERNEL_TYPE>
{
public:
	typedef typename SparseMatrix<T,KERNEL_TYPE>::ROW ROW;
	struct SPOP
	{
	public:
		SPOP(const FixedSparseMatrix<T,KERNEL_TYPE>& matrix):_matrix(&matrix),_eigenMatrix(NULL){}
		SPOP(const Eigen::SparseMatrix<T,0,sizeType>& matrix):_matrix(NULL),_eigenMatrix(&matrix){}
		void operator()(const ROW& rhs,ROW& result) const
		{
			//put
			if(_matrix)putMatrix(rhs,result);
			else putEigenMatrix(rhs,result);
			//compact
			if(result.empty())return;
			std::sort(result.begin(),result.end(),typename SparseMatrix<T,KERNEL_TYPE>::SortByCol());
			sizeType nz=0;
			for(sizeType i=1;i<(sizeType)result.size();i++)
				if(result[i].second == result[nz].second)
					result[nz].first+=result[i].first;
				else result[++nz]=result[i];
			result.resize(++nz);
		}
		void putMatrix(const ROW& rhs,ROW& result) const
		{
			const vector<sizeType>& rOff=_matrix->getRowOffset();
			const ROW& rVal=_matrix->getValue();
			result.clear();
			for(sizeType i=0;i<(sizeType)rhs.size();i++)
			for(sizeType beg=rOff[rhs[i].second],end=rOff[rhs[i].second+1];beg!=end;beg++)
				result.push_back(std::pair<T,sizeType>(rVal[beg].first*rhs[i].first,rVal[beg].second));
		}
		void putEigenMatrix(const ROW& rhs,ROW& result) const
		{
			result.clear();
			for(sizeType i=0;i<(sizeType)rhs.size();i++)
			for(typename Eigen::SparseMatrix<T,0,sizeType>::InnerIterator it(*_eigenMatrix,rhs[i].second);it;++it)
				result.push_back(std::pair<T,sizeType>(it.value()*rhs[i].first,it.row()));
		}
		sizeType n() const{return _matrix ? _matrix->rows() : _eigenMatrix->rows();}
		const FixedSparseMatrix<T,KERNEL_TYPE>* _matrix;
		const Eigen::SparseMatrix<T,0,sizeType>* _eigenMatrix;
	};
public:
	typedef typename KERNEL_TYPE::Vec Vec;
    using Solver<T,KERNEL_TYPE>::_residualOut;
    using Solver<T,KERNEL_TYPE>::_iterationsOut;
	SymAINVPreconSolver(){setSolverParameter(0.01f,0.01f);}
	virtual void setSolverParameter(T tol1,T tol2){_tol1=tol1;_tol2=tol2;}
    virtual void setMatrix(const FixedSparseMatrix<T,KERNEL_TYPE>& matrix,bool syncPrecon){computeAInv(SPOP(matrix));}
    virtual void setMatrix(const Eigen::SparseMatrix<T,0,sizeType>& matrix,bool syncPrecon){computeAInv(SPOP(matrix));}
    typename Solver<T, KERNEL_TYPE>::SOLVER_RESULT solve(const Vec& rhs,Vec& result)
	{
		_ZInvDZT.multiply(rhs,result);
		return Solver<T, KERNEL_TYPE>::SUCCESSFUL;
	}
private:
	template <typename SPOP_TYPE>
	void computeAInv(SPOP_TYPE op)
	{
		//init memory and set Z=I
		ROW V,AX,tmp;
		SparseMatrix<T,KERNEL_TYPE> ZTTmp(op.n());
		vector<boost::shared_ptr<ROW> >& rows=ZTTmp.getValue();
		for(sizeType i=0;i<op.n();i++)
			ZTTmp.addToElement(i,i,1.0f);
		Vec invD(op.n());

		//main loop
		for(sizeType i=0;i<op.n();i++)
		{
			if(i%100 == 0)
			INFOV("SymAINV computed row: %d!",i)

			//Step 1: V=drop(AZ)
			drop(*(rows[i]),_tol1);
			op(*(rows[i]),V);

			//Step 2: p_j^{i-1}=<v,z_j^{i-1}>
			//Step 3: z_j-=drop((p_j/p_i)z_i)
			invD[i]=1.0f/dot(*(ZTTmp.getValue()[i]),V);
			ASSERT_MSG(dot(*(ZTTmp.getValue()[i]),V) > 1E-10f,"AINV BreakDown!")
			//case from i+1
			OMP_PARALLEL_FOR_I(OMP_PRI(AX,tmp))
			for(sizeType j=i+1;j<op.n();j++)
			{
				T coef=dot(*(ZTTmp.getValue()[j]),V)*invD[i];
				AX=*(rows[i]);
				scale(AX,-coef);
				drop(AX,_tol2);
				add(AX,*(rows[j]),tmp);
			}
		}

		//assemble final Z
		FixedSparseMatrix<T,KERNEL_TYPE> ZT,Z;
		ZT.constructFromMatrix(ZTTmp);
		ZT.multiplyRow(invD.cwiseSqrt());
		ZT.transpose(Z);
		Z.multiplyMatrix(ZT,_ZInvDZT);
	}
	static void drop(ROW& A,T tol)
	{
		if(tol == 0.0f)return;
		sizeType j=0;
		for(sizeType i=0;i<(sizeType)A.size();i++)
		if(fabs(A[i].first) > tol)A[j++]=A[i];
		A.resize(j);
	}
	static T dot(const ROW& A,const ROW& B)
	{
		sizeType i=0,j=0;
		T out=0.0f;
		while(i < (sizeType)A.size() && j < (sizeType)B.size())
		{
			if(A[i].second < B[j].second)i++;
			else if(A[i].second > B[j].second)j++;
			else{out+=A[i].first*B[j].first;i++;j++;}
		}
		return out;
	}
	static void scale(ROW& A,T coef)
	{
		for(sizeType i=0;i<(sizeType)A.size();i++)
			A[i].first*=coef;
	}
	static void add(const ROW& A,ROW& B,ROW& tmp)
	{
		tmp=B;
		B.clear();
		sizeType i=0,j=0;
		while(i < (sizeType)A.size() && j < (sizeType)tmp.size())
		{
			if(A[i].second < tmp[j].second){B.push_back(A[i]);i++;}
			else if(A[i].second > tmp[j].second){B.push_back(tmp[j]);j++;}
			else{B.push_back(std::pair<T,sizeType>(A[i].first+tmp[j].first,A[i].second));i++;j++;}
		}
		while(i < (sizeType)A.size()){B.push_back(A[i]);i++;}
		while(j < (sizeType)tmp.size()){B.push_back(tmp[j]);j++;}
	}
	T _tol1,_tol2;
	FixedSparseMatrix<T,KERNEL_TYPE> _ZInvDZT;
};

PRJ_END

#endif
