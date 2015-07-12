#ifndef MAT_VEC_H
#define MAT_VEC_H

#include "Config.h"
#include "MathBasic.h"
#include "IO.h"

#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>

PRJ_BEGIN

//============================================================================
// Vector Interface
template <typename T>
struct Kernel {
    template<typename T2>
    struct Rebind {
        typedef Kernel<T2> value;
    };
	typedef T Scalar;
    typedef Eigen::Matrix<T,-1,1> Vec;
	typedef Eigen::Matrix<T,-1,-1> Mat;
    static FORCE_INLINE T dot(const Vec& x, const Vec& y,sizeType n=-1) {
		if(n == -1)return x.dot(y);
		else return x.block(0,0,n,1).dot(y.block(0,0,n,1));
	}
    static FORCE_INLINE T norm(const Vec& x,sizeType n=-1) {
        if(n == -1)return x.norm();
		else return x.block(0,0,n,1).norm();
    }
    static FORCE_INLINE T norm1(const Vec& x,sizeType n=-1) {
        if(n == -1)return x.cwiseAbs().sum();
		else return x.block(0,0,n,1).cwiseAbs().sum();
    }
    static FORCE_INLINE void copy(const Vec& x,Vec& y,sizeType n=-1) {
        if(n == -1)y=x;
		else y.block(0,0,n,1)=x.block(0,0,n,1);
	}
    static FORCE_INLINE void copy(const Mat& x,Mat& y) {
        y=x;
	}
    static FORCE_INLINE void ncopy(const Vec& x,Vec& y,sizeType n=-1) {
		if(n == -1)y=-x;
		else y.block(0,0,n,1)=-x.block(0,0,n,1);
    }
    static FORCE_INLINE sizeType indexAbsMax(const Vec& x,sizeType n=-1) {
        if(n == -1)n=x.size();
		sizeType maxInd = 0;
        T maxValue = (T)0.0f;
        for(sizeType i = 0; i < n; ++i) {
            if(std::abs(x[i]) > maxValue) {
                maxValue = std::abs(x[i]);
                maxInd = i;
            }
        }
        return maxInd;
    }
    static FORCE_INLINE T absMax(const Vec& x,sizeType n=-1) {
        if(n == -1)return std::abs(x[indexAbsMax(x)]);
		else return std::abs(x.block(0,0,n,1)[indexAbsMax(x.block(0,0,n,1))]);
	}
    static FORCE_INLINE void scale(T alpha,Vec& y,sizeType n=-1) {
		if(n == -1)y*=alpha;
		else y.block(0,0,n,1)*=alpha;
	}
    static FORCE_INLINE void addScaled(const T& alpha, const Vec& x, Vec& y,sizeType n=-1) {
		if(n == -1)y+=x*alpha;
		else y.block(0,0,n,1)+=x.block(0,0,n,1)*alpha;
	}
    static FORCE_INLINE void add(const Vec& x,const Vec& y,Vec& result,sizeType n=-1) {
        if(n == -1)result=x+y;
		else result.block(0,0,n,1)=x.block(0,0,n,1)+y.block(0,0,n,1);
	}
    static FORCE_INLINE void sub(const Vec& x,const Vec& y,Vec& result,sizeType n=-1) {
        if(n == -1)result=x-y;
		else result.block(0,0,n,1)=x.block(0,0,n,1)-y.block(0,0,n,1);
    }
    static FORCE_INLINE void cwiseProd(const Vec& x, const Vec& y, Vec& z,sizeType n=-1) {
		if(n == -1)z.array()=x.array()*y.array();
		else z.block(0,0,n,1).array()=
			x.block(0,0,n,1).array()*
			y.block(0,0,n,1).array();
	}
    static FORCE_INLINE void cwisePow(const T& p, Vec& x,sizeType n=-1) {
		if(n == -1)x=x.cwisePow(p);
		else x.block(0,0,n,1)=x.block(0,0,n,1).cwisePow(p);
	}
    static FORCE_INLINE void zero(Vec& v,sizeType n=-1) {set(v,0.0f,n);}
    static FORCE_INLINE void set(Vec& v,const T& val,sizeType n=-1) {
        if(n == -1)v.setConstant(val);
		else v.block(0,0,n,1).setConstant(val);
    }
	static FORCE_INLINE void mul(const Mat& A,const Vec& x,Vec& out){
		out=A*x;
	}
	static FORCE_INLINE void mulT(const Mat& A,const Vec& x,Vec& out){
		out=A.transpose()*x;
	}
	static FORCE_INLINE sizeType rows(const Mat& A){
		return A.rows();
	}
	static FORCE_INLINE sizeType cols(const Mat& A){
		return A.cols();
	}
	static FORCE_INLINE T colNorm(const Mat& A,sizeType i){
		return A.col(i).norm();
	}
	template <typename A,typename B>
	static FORCE_INLINE void addMulT(A& x,const T& b,const B& c){x+=b*c;}
	template <typename A,typename B>
	static FORCE_INLINE void subMulT(A& x,const T& b,const B& c){x-=b*c;}
	static FORCE_INLINE T absV(const T& val){return std::abs(val);}
	static FORCE_INLINE T ZeroV(){return (T)0.0f;}
};

template <typename T,int R,int C>
struct Kernel<Eigen::Matrix<T,R,C> > {
	typedef T Scalar;
    typedef typename Eigen::Matrix<T,-1,1> Vec;
	template <typename A,typename B>
	static FORCE_INLINE void addMulT(A& x,const typename Eigen::Matrix<T,R,C>& b,const B& c){x+=b.transpose()*c;}
	template <typename A,typename B>
	static FORCE_INLINE void subMulT(A& x,const typename Eigen::Matrix<T,R,C>& b,const B& c){x-=b.transpose()*c;}
	static FORCE_INLINE T absV(const Eigen::Matrix<T,R,C>& val){return val.norm();}
	static FORCE_INLINE Eigen::Matrix<T,R,C> ZeroV(){return Eigen::Matrix<T,R,C>::Zero();}
};

//============================================================================
// Sparse Matrix Iterator
template <typename T>
struct ConstSMIterator {
public:
    ConstSMIterator():_r(-1) {}
    ConstSMIterator(sizeType r,typename std::vector<std::pair<T,sizeType> >::const_iterator it):_iter(it),_r(r) {}
    virtual ~ConstSMIterator() {}
    ConstSMIterator& operator++() {
        _iter++;
        return *this;
    }
    bool operator==(const ConstSMIterator& other) const {
        if(_r == other._r) {
            if(_r == -1)
                return true;
            return _iter == other._iter;
        }
        return false;
    }
    bool operator!=(const ConstSMIterator& other) const {
        return !operator==(other);
    }
    const T& operator*() const {
        return _iter->first;
    }
    sizeType col() const {
        return _iter->second;
    }
    sizeType row() const {
        return _r;
    }
    std::pair<T,sizeType> pair() const {
        return *_iter;
    }
protected:
    typename std::vector<std::pair<T,sizeType> >::const_iterator _iter;
    sizeType _r;
};
template <typename T>
struct SMIterator : public ConstSMIterator<T> {
public:
    using ConstSMIterator<T>::_iter;
    SMIterator() {}
    SMIterator(sizeType r,typename std::vector<std::pair<T,sizeType> >::const_iterator it):ConstSMIterator<T>(r,it) {}
    SMIterator& operator++() {
        _iter++;
        return *this;
    }
    T& operator*() {
        return (T&)(static_cast<ConstSMIterator<T>&>(*this).operator*());
    }
};

//============================================================================
// CSC Sparse Matrix
template <typename T,typename KERNEL_TYPE>
struct SparseMatrix : public Serializable{
public:
    typedef typename std::pair<T,sizeType> ROWE;
    typedef typename std::vector<ROWE> ROW;
    struct SortByCol {
        bool operator()(const ROWE& A,const ROWE& B) const {
            return A.second < B.second;
        }
    };
    typedef typename KERNEL_TYPE::Vec Vec;
    SparseMatrix(sizeType n=0, sizeType nz=7, sizeType m=-1) 
	:Serializable(-1){reset(n,nz,m);}
    void reset(sizeType n=0, sizeType nz=7, sizeType m=-1) {
        if(m==-1)m=n;

        resize(n,m);
        for(sizeType i=0; i<_n; ++i) {
            _value[i].reset(new ROW());
            _value[i]->reserve(nz);
        }
    }
    void clear(void) {
        _n=0;
        _m=0;
        _value.clear();
    }
    void zero(void) {
        for(sizeType i=0; i<_n; ++i)
            _value[i].reset((ROW*)0);
    }
    void resize(sizeType n,sizeType m=-1) {
        _n=n;
        if(m==-1)m=n;
        _m=m;
        _value.resize(n);
        zero();
    }
    sizeType rows() const {
        return _n;
    }
    sizeType cols() const {
        return _m;
    }
    sizeType numElement(const sizeType& i) const {
        if(!(_value[i]))
            return 0;
        else return _value[i]->size();
    }
    T operator()(sizeType i, sizeType j) const {
        ASSERT(j < _m)
        if(!(_value[i]))
            return (T)0.0f;
        const ROW& r=*(_value[i]);
        for(sizeType k=0; k<(sizeType)r.size(); ++k) {
            if(r[k].second==j)
                return r[k].first;
            else if(r[k].second>j)
                return (T)0.0f;
        }
        return (T)0.0f;
    }
    static bool rowEquals(ConstSMIterator<T> beg,ConstSMIterator<T> end,
                          ConstSMIterator<T> oBeg,ConstSMIterator<T> oEnd,
                          const T& eps=ScalarUtil<T>::scalar_eps) {
        while(beg != end && oBeg != oEnd) {
            if(oBeg.col() < beg.col()) {
                if(std::abs(*oBeg) > eps)return false;
                ++oBeg;
            } else if(beg.col() < oBeg.col()) {
                if(std::abs(*beg) > eps)return false;
                ++beg;
            } else {
                if(std::abs(*beg-*oBeg) > eps)return false;
                ++beg;
                ++oBeg;
            }
        }
        while(beg != end) {
            if(std::abs(*beg) > eps)return false;
            ++beg;
        }
        while(oBeg != oEnd) {
            if(std::abs(*oBeg) > eps)return false;
            ++oBeg;
        }
        return true;
    }
    bool equals(const SparseMatrix& other,const T& eps=ScalarUtil<T>::scalar_eps) const {
        if(_n != other._n || _m != other._m)
            return false;
        for(sizeType i=0; i<_n; i++) {
            ConstSMIterator<T> cBeg=begin(i);
            ConstSMIterator<T> cEnd=end(i);
            ConstSMIterator<T> oBeg=other.begin(i);
            ConstSMIterator<T> oEnd=other.end(i);
            if(!rowEquals(cBeg,cEnd,oBeg,oEnd,eps))
                return false;
        }
        return true;
    }
    bool isZero(const T& eps=ScalarUtil<T>::scalar_eps) const {
        for(sizeType i=0; i<_n; i++)
            if((_value[i])) {
                const ROW& r=*(_value[i]);
                const sizeType nrE=(sizeType)r.size();
                for(sizeType j=0; j<nrE; j++)
                    if(std::abs(r[j].first) > eps)
                        return false;
            }
        return true;
    }
    void setElement(sizeType i,sizeType j,T newValue) {
        ASSERT(j < _m)
        if(!(_value[i]))
            _value[i].reset(new ROW());
        ROW& r=*(_value[i]);
        for(sizeType k=0; k<(sizeType)r.size(); ++k) {
            if(r[k].second==j) {
                r[k].first=newValue;
                return;
            } else if(r[k].second>j) {
                r.insert(r.begin()+k,std::make_pair(newValue,j));
                return;
            }
        }
        r.push_back(std::make_pair(newValue,j));
    }
    void addToElement(sizeType i,sizeType j,T incValue) {
        ASSERT(j < _m)
        sizeType k=0;
        if(!(_value[i]))
            _value[i].reset(new ROW());
        ROW& r=*(_value[i]);
        for(; k<(sizeType)r.size(); ++k) {
            if(r[k].second==j) {
                r[k].first+=incValue;
                return;
            } else if(r[k].second>j) {
                r.insert(r.begin()+k,std::make_pair(incValue,j));
                return;
            }
        }
        _value[i]->push_back(std::make_pair(incValue,j));
    }
    void add(const SparseMatrix& other,const T& coef=1.0f,const T& eps=ScalarUtil<T>::scalar_eps) {
        ASSERT(other._n == _n && other._m == _m)
        for(sizeType i=0; i<_n; i++)
            if((other._value[i]))
                addSparseRowC(i,*(other._value[i]),coef,eps);
    }
    void sub(const SparseMatrix& other,const T& eps=ScalarUtil<T>::scalar_eps) {
        add(other,-1.0f,eps);
    }
    void addSparseRow(sizeType i,const ROW &r,const T& eps=ScalarUtil<T>::scalar_eps) {
        if(!(_value[i])) {
            _value[i].reset(new ROW());
            *(_value[i])=r;
            return;
        }
        boost::shared_ptr<ROW> newR(new ROW);
        {
            const ROW& ro=*(_value[i]);
            typename ROW::const_iterator rBeg=r.begin();
            typename ROW::const_iterator rEnd=r.end();
            typename ROW::const_iterator roBeg=ro.begin();
            typename ROW::const_iterator roEnd=ro.end();
            merge(rBeg,rEnd,roBeg,roEnd,*newR);
        }
        if(newR->empty())
            newR.reset((ROW*)0);
        _value[i]=newR;
    }
    void addSparseRow(sizeType i,const std::vector<sizeType>& index,const std::vector<T>& value,const T& eps=ScalarUtil<T>::scalar_eps) {
        ROW tmp;
        for(sizeType j=0; j<(sizeType)index.size(); j++)
            tmp.push_back(std::make_pair(value[j],index[j]));
        addSparseRow(i,tmp,eps);
    }
    void addSparseRowC(sizeType i,const ROW& r,const T& coef,const T& eps=ScalarUtil<T>::scalar_eps) {
        ROW valuesM=r;
        for(typename ROW::iterator beg=valuesM.begin(),end=valuesM.end(); beg!=end; beg++)
            (*beg).first*=coef;
        addSparseRow(i,valuesM,eps);
    }
    void addSparseRowC(sizeType i,const std::vector<sizeType>& index,const std::vector<T>& value,const T& coef,const T& eps=ScalarUtil<T>::scalar_eps) {
        std::vector<T> valuesM=value;
        for(typename std::vector<T>::iterator beg=valuesM.begin(),end=valuesM.end(); beg!=end; beg++)
            (*beg)*=coef;
        addSparseRow(i,index,valuesM,eps);
    }
    void addSparseRowC(sizeType i,ConstSMIterator<T> beg,ConstSMIterator<T> end,const T& coef,const T& eps=ScalarUtil<T>::scalar_eps) {
        ROW valuesM;
        for(; beg!=end; ++beg)
            valuesM.push_back(std::make_pair((*beg)*coef,beg.col()));
        addSparseRow(i,valuesM,eps);
    }
    static void merge(typename ROW::const_iterator rBeg,
                      typename ROW::const_iterator rEnd,
                      typename ROW::const_iterator roBeg,
                      typename ROW::const_iterator roEnd,ROW& newR,const T& eps=ScalarUtil<T>::scalar_eps) {
        while(rBeg != rEnd && roBeg != roEnd) {
            const ROWE& e=*rBeg;
            const ROWE& eo=*roBeg;
            if(e.second < eo.second) {
                newR.push_back(e);
                ++rBeg;
            } else if(e.second > eo.second) {
                newR.push_back(eo);
                ++roBeg;
            } else {
                T tmp=e.first+eo.first;
                if(std::abs(tmp) > eps)
                    newR.push_back(std::make_pair(tmp,e.second));
                ++rBeg;
                ++roBeg;
            }
        }
        while(rBeg != rEnd) {
            newR.push_back(*rBeg);
            ++rBeg;
        }
        while(roBeg != roEnd) {
            newR.push_back(*roBeg);
            ++roBeg;
        }
    }
    void symmetricRemoveRolCol(sizeType i) {
        if((_value[i])) {
            ROW& r=*(_value[i]);
            const sizeType nrE=(sizeType)r.size();
            for(sizeType a=0; a<(sizeType)nrE; ++a) {
                sizeType j=r[a].second;
                if(!(_value[j]))
                    continue;

                ROW& ro=*(_value[j]);
                const sizeType nrEO=(sizeType)ro.size();
                for(sizeType b=0; b<(sizeType)nrEO; ++b) {
                    if(ro[b].second==i) {
                        ro.erase(ro.begin()+b);
                        break;
                    }
                }
            }
        }
        _value[i].reset((ROW*)0);
    }
    bool isSymmetric(const T& tol=ScalarUtil<T>::scalar_eps) const {
        for(sizeType i=0; i<_n; i++) {
            if(!(_value[i]))
                continue;
            const ROW& r=*(_value[i]);
            const sizeType nrE=(sizeType)r.size();
            for(sizeType j=0; j<nrE; j++)
                if(std::abs(operator()(i,r[j].second)-operator()(r[j].second,i)) > tol)
                    return false;
        }
        return true;
    }
    bool isDiagonalDominant() const {
        for(sizeType i=0; i<_n; i++) {
            if(!(_value[i]))
                continue;

            T res=(T)0.0f;
            const ROW& r=*(_value[i]);
            const sizeType nrE=(sizeType)r.size();
            for(sizeType j=0; j<(sizeType)nrE; j++) {
                if(r[j].second == i)
                    res+=std::abs(r[j].first);
                else res-=std::abs(r[j].first);
            }
            if(res < (T)0.0f)
                return false;
        }
        return true;
    }
    //io
    void toEigen(const sizeType& r0,const sizeType& c0,std::vector<Eigen::Triplet<T,sizeType> >& m,bool sym=false) const {
        for(sizeType i=0; i<_n; ++i) {
            if(!(_value[i]))
                continue;
            const ROW& r=*(_value[i]);
            const sizeType nrE=(sizeType)r.size();
            for(sizeType j=0; j<nrE; ++j) {
                m.push_back(Eigen::Triplet<T,sizeType>((sizeType)(r0+i),(sizeType)(c0+r[j].second),r[j].first));
                if(sym)
                    m.push_back(Eigen::Triplet<T,sizeType>((sizeType)(c0+r[j].second),(sizeType)(r0+i),r[j].first));
            }
        }
    }
    void toEigen(Eigen::SparseMatrix<T,0,sizeType>& m) const {
        m.resize(_n,_m);
        std::vector<Eigen::Triplet<T,sizeType> > trips;
        toEigen(0,0,trips);
        m.setFromTriplets(trips.begin(),trips.end());
    }
    void fromEigen(const Eigen::SparseMatrix<T,0,sizeType>& m) {
        resize(m.rows(),m.cols());
        zero();
        for(sizeType k=0; k<m.outerSize(); ++k)
            for(typename Eigen::SparseMatrix<T,0,sizeType>::InnerIterator it(m,k); it; ++it)
                setElement(it.row(),it.col(),it.value());
    }
    bool write(std::ostream& output) const {
        writeBinaryData(_n,output);
        writeBinaryData(_m,output);
        for(sizeType i=0; i<_n; i++) {
            sizeType nr=0;
            if((_value[i]))
                nr=(sizeType)_value[i]->size();
            writeBinaryData(nr,output);
            if(nr > 0)
                writeVector(*(_value[i]),output);
        }
		return output.good();
    }
    bool read(std::istream& input) {
        readBinaryData(_n,input);
        readBinaryData(_m,input);
        _value.resize(_n);
        for(sizeType i=0; i<_n; i++) {
            sizeType nr;
            readBinaryData(nr,input);
            if(nr == 0)
                _value[i].reset((ROW*)0);
            else {
                _value[i].reset(new ROW());
                readVector(*(_value[i]),input);
            }
        }
		return input.good();
    }
    //mat vec
    template <typename VEC,typename VEC_OUT>
    void multiply(const VEC& x,VEC_OUT& result) const {
        //ASSERT(_m==x.size() && _n==result.size());
        //result.resize(_n);
        for(sizeType i=0; i<_n; ++i) {
            if(!(_value[i]))
                continue;
            const ROW& r=*(_value[i]);
            const sizeType nrE=(sizeType)r.size();

            result[i]=0;
            for(sizeType j=0; j<nrE; ++j)
                result[i]+=r[j].first*x[r[j].second];
        }
    }
    template <typename VEC,typename VEC_OUT>
    void multiplySubtract(const VEC& x,VEC_OUT& result) const {
        //ASSERT(_m==x.size() && _n==result.size());
        //result.resize(_n);
        for(sizeType i=0; i<_n; ++i) {
            if(!(_value[i]))
                continue;
            const ROW& r=*(_value[i]);
            const sizeType nrE=(sizeType)r.size();

            for(sizeType j=0; j<nrE; ++j)
                result[i]-=r[j].first*x[r[j].second];
        }
    }
    template <typename VEC,typename VEC_OUT>
    void multiplyAdd(const VEC& x,VEC_OUT& result) const {
        //ASSERT(_m==x.size() && _n==result.size());
        //result.resize(_n);
        for(sizeType i=0; i<_n; ++i) {
            if(!(_value[i]))
                continue;
            const ROW& r=*(_value[i]);
            const sizeType nrE=(sizeType)r.size();

            for(sizeType j=0; j<nrE; ++j)
                result[i]+=r[j].first*x[r[j].second];
        }
    }
    //mat vec transposed
    template <typename VEC,typename VEC_OUT>
    void multiplyTranspose(const VEC& x,VEC_OUT& result) const {
        //ASSERT(_m==x.size() && _n==result.size());
        //result.resize(_n);
        OMP_PARALLEL_FOR_
        for(sizeType i=0; i<_m; i++)
            result[i]=(T)0.0f;
        multiplyTransposeAdd(x,result);
    }
    template <typename VEC,typename VEC_OUT>
    void multiplyTransposeSubtract(const VEC& x,VEC_OUT& result) const {
        //ASSERT(_m==x.size() && _n==result.size());
        //result.resize(_n);
        for(sizeType i=0; i<_n; ++i) {
            if(!(_value[i]))
                continue;
            const ROW& r=*(_value[i]);
            const sizeType nrE=(sizeType)r.size();

            for(sizeType j=0; j<nrE; ++j)
                result[r[j].second]-=r[j].first*x[i];
        }
    }
    template <typename VEC,typename VEC_OUT>
    void multiplyTransposeAdd(const VEC& x,VEC_OUT& result) const {
        //ASSERT(_m==x.size() && _n==result.size());
        //result.resize(_n);
        for(sizeType i=0; i<_n; ++i) {
            if(!(_value[i]))
                continue;
            const ROW& r=*(_value[i]);
            const sizeType nrE=(sizeType)r.size();

            for(sizeType j=0; j<nrE; ++j)
                result[r[j].second]+=r[j].first*x[i];
        }
    }
    //mat mat
    void multiplyMatrix(const SparseMatrix& b,SparseMatrix& out) const {
        out.reset(_n,0,b._m);

        sizeType nnz=0;
        std::vector<sizeType> off(b._m,-1);
        std::vector<sizeType> col(b._m,-1);
        for(sizeType i=0; i<_n; i++) {
            boost::shared_ptr<ROW> oPtr(new ROW());
            ROW& rOut=*(oPtr);

            if(!(_value[i]))
                continue;
            const ROW& r=*(_value[i]);
            const sizeType nrE=(sizeType)r.size();
            for(sizeType j=0; j<nrE; j++) {
                ROWE c=r[j];
                if(!(b._value[c.second]))
                    continue;
                const ROW& rB=*(b._value[c.second]);
                const sizeType nrEB=(sizeType)rB.size();
                for(sizeType k=0; k<nrEB; k++) {
                    ROWE cB=rB[k];
                    ROWE cOut=std::make_pair(c.first*cB.first,cB.second);
                    if(off[cOut.second] == -1) {
                        off[cOut.second]=nnz;
                        col[nnz]=cOut.second;
                        rOut.push_back(cOut);
                        nnz++;
                    } else {
                        rOut[off[cOut.second]].first+=cOut.first;
                    }
                }
            }

            if(!rOut.empty())
                out._value[i]=oPtr;
            for(sizeType j=0; j<nnz; j++)
                off[col[j]]=-1;
            nnz=0;
            std::sort(oPtr->begin(),oPtr->end(),SparseMatrix<T,KERNEL_TYPE>::SortByCol());
        }
    }
    void multiplyMatrixTransposed(const SparseMatrix& b,SparseMatrix& out) const {
        SparseMatrix bTmp;
        b.transpose(bTmp);
        multiplyMatrix(bTmp,out);
    }
    void transposedMultiplyMatrix(const SparseMatrix& b,SparseMatrix& out) const {
        SparseMatrix ATmp;
        transpose(ATmp);
        ATmp.multiplyMatrix(b,out);
    }
    //transpace
    void transpose(SparseMatrix& b) const {
        b.reset(_m,0,_n);
        for(sizeType i=0; i<_n; i++) {
            if(!(_value[i]))
                continue;
            const ROW& r=*(_value[i]);
            const sizeType nrE=(sizeType)r.size();

            for(sizeType j=0; j<nrE; j++)
                b.setElement(r[j].second,i,r[j].first);
        }
    }
    //iterator
    ConstSMIterator<T> begin(sizeType r) const {
        if((_value[r]))
            return ConstSMIterator<T>(r,_value[r]->begin());
        else return ConstSMIterator<T>(-1,typename ROW::iterator());
    }
    ConstSMIterator<T> end(sizeType r) const {
        if((_value[r]))
            return ConstSMIterator<T>(r,_value[r]->end());
        else return ConstSMIterator<T>(-1,typename ROW::iterator());
    }
    SMIterator<T> begin(sizeType r) {
        if((_value[r]))
            return SMIterator<T>(r,_value[r]->begin());
        else return SMIterator<T>(-1,typename ROW::iterator());
    }
    SMIterator<T> end(sizeType r) {
        if((_value[r]))
            return SMIterator<T>(r,_value[r]->end());
        else return SMIterator<T>(-1,typename ROW::iterator());
    }
    FORCE_INLINE const std::vector< boost::shared_ptr<ROW> >& getValue() const{return _value;}
    FORCE_INLINE std::vector< boost::shared_ptr<ROW> >& getValue(){return _value;}
private:
    sizeType _n;
    sizeType _m;
    std::vector< boost::shared_ptr<ROW> > _value;
};

template <typename T,typename KERNEL_TYPE>
struct FixedSparseMatrix : public Serializable{
public:
	typedef typename KERNEL_TYPE::Scalar Scalar;
    typedef typename std::pair<T,sizeType> ROWE;
    typedef typename std::vector<ROWE> ROW;
    typedef typename KERNEL_TYPE::Vec Vec;
    struct Lss
    {
        bool operator()(const Eigen::Triplet<T,sizeType>& A,const Eigen::Triplet<T,sizeType>& B) const
        {return A.row() < B.row() || (A.row() == B.row() && A.col() < B.col());}
    };
    FixedSparseMatrix(sizeType n=0,sizeType m=-1)
	:Serializable(-1),_n(n),_m(m==-1 ? n : m),_value(0),_rowStart(n+1) {}
	virtual ~FixedSparseMatrix(){}
    void clear(void) {
        _n=0;
        _m=0;
        _value.clear();
        _rowStart.clear();
    }
    void zero(void) {
        _value.clear();
        _rowStart.assign(_n+1,0);
    }
    void reset(sizeType n,sizeType nnz,sizeType m) {
        resize(n,m);
        _value.reserve(nnz*_n);
    }
    void resize(sizeType n,sizeType m=-1) {
        _n=n;
        _m=m==-1 ? n : m;
        _rowStart.resize(n+1);
        _value.clear();
        zero();
    }
    template <typename ITER>
    void buildFromTriplets(ITER beg,ITER end) {
        sizeType nnz=0;
        _rowStart.assign(_n+1,0);
        {
            ITER begT=beg;
            ITER endT=end;
            for(; begT!=endT; begT++)
            {
                if(begT->row() < 0)
                    continue;
                _rowStart[begT->row()+1]++;
                nnz++;
            }
            for(sizeType i=0; i<_n; i++)
                _rowStart[i+1]+=_rowStart[i];
        }
        {
            std::vector<sizeType> tmp=_rowStart;
            _value.resize(nnz);
            ITER begT=beg;
            ITER endT=end;
            for(; begT!=endT; begT++) 
            {
                if(begT->row() < 0)
                    continue;
                sizeType& r=tmp[begT->row()];
                _value[r]=std::make_pair(begT->value(),begT->col());
                r++;
            }
        }
        {
            for(sizeType i=0; i<_n; i++)
                std::sort(_value.begin()+_rowStart[i],
                          _value.begin()+_rowStart[i+1],
                          typename SparseMatrix<T,KERNEL_TYPE>::SortByCol());
        }
    }
    void buildFromTripletsDepulicate(std::vector<Eigen::Triplet<T,sizeType> >& trips,Scalar eps=1E-8f)
    {
        std::sort(trips.begin(),trips.end(),Lss());
        
        const sizeType nr=(sizeType)trips.size();
        _rowStart.assign(_n+1,0);
        _value.reserve(nr);
        _value.clear();

        sizeType r=-1,c=-1;
        T curr=KERNEL_TYPE::ZeroV();
        for(sizeType i=0;i<nr;i++)
        {
            const Eigen::Triplet<T,sizeType>& trip=trips[i];
            if(r == trip.row() && c == trip.col())
                curr+=trip.value();
            else{
                if(r >= 0 && c >= 0 && (eps == (Scalar)0.0f || KERNEL_TYPE::absV(curr) > eps)){
                    //feed value
                    _rowStart[r+1]++;
                    _value.push_back(make_pair(curr,c));
                }
                //new entry initialize
                r=trip.row();
                c=trip.col();
                curr=trip.value();
            }
        }
        if(r >= 0 && c >= 0 && KERNEL_TYPE::absV(curr) > eps){
            //feed value
            _rowStart[r+1]++;
            _value.push_back(make_pair(curr,c));
        }
        //build row offset
        for(sizeType i=1;i<=_n;i++)
            _rowStart[i]+=_rowStart[i-1];
    }
    sizeType rows() const {
        return _n;
    }
    sizeType cols() const {
        return _m;
    }
    sizeType numElement(const sizeType& i) const {
        return _rowStart[i+1]-_rowStart[i];
    }
    void constructFromMatrix(const SparseMatrix<T,KERNEL_TYPE>& matrix) {
        resize(matrix.rows(),matrix.cols());
        _rowStart[0]=0;
        for(sizeType i=0; i<_n; ++i) {
            sizeType sz=matrix.numElement(i);
            _rowStart[i+1]=(sizeType)(_rowStart[i]+sz);
        }

        _value.resize(_rowStart[_n]);
        typename ROW::iterator vIter=_value.begin();
        for(sizeType i=0; i<_n; ++i) {
            ConstSMIterator<T> beg=matrix.begin(i),end=matrix.end(i);
            for(; beg!=end; ++beg) {
                *vIter=beg.pair();
                vIter++;
            }
        }
    }
    void toSparseMatrix(SparseMatrix<T,KERNEL_TYPE>& matrix) {
        matrix.resize(_n,_m);
        for(sizeType i=0; i<_n; i++) {
            ROW tmp;
            for(sizeType j=_rowStart[i]; j<_rowStart[i+1]; j++)
                tmp.push_back(_value[j]);
            matrix.addSparseRow(i,tmp);
        }
    }
    //io
	template <typename ST>
    void toEigen(const sizeType& r,const sizeType& c,std::vector<Eigen::Triplet<T,ST> >& m,bool sym=false) const {
        for(sizeType i=0; i<_n; ++i) {
            for(sizeType j=_rowStart[i]; j<_rowStart[i+1]; ++j) {
                m.push_back(Eigen::Triplet<T,ST>((ST)(r+i),(ST)(c+_value[j].second),_value[j].first));
                if(sym)
                    m.push_back(Eigen::Triplet<T,ST>((ST)(c+_value[j].second),(ST)(r+i),_value[j].first));
            }
        }
    }
    template <int O,typename ST>
	void toEigen(Eigen::SparseMatrix<T,O,ST>& m) const {
        m.resize((ST)_n,(ST)_m);
        std::vector<Eigen::Triplet<T,ST> > trips;
        toEigen<ST>(0,0,trips);
        m.setFromTriplets(trips.begin(),trips.end());
    }
	template <int O,typename ST>
    void fromEigen(const Eigen::SparseMatrix<T,O,ST>& m) {
        clear();
        resize(m.rows(),m.cols());
		std::vector<Eigen::Triplet<T,sizeType> > trips;
        for(sizeType k=0; k<m.outerSize(); ++k)
            for(typename Eigen::SparseMatrix<T,O,ST>::InnerIterator it(m,(ST)k); it; ++it)
                trips.push_back(Eigen::Triplet<T,sizeType>(it.row(),it.col(),it.value()));
		buildFromTriplets(trips.begin(),trips.end());
	}
    bool write(std::ostream& output) const {
        writeBinaryData(_n,output);
        writeBinaryData(_m,output);
        writeVector(_value,output);
        writeVector(_rowStart,output);
		return output.good();
    }
    bool read(std::istream& input) {
        readBinaryData(_n,input);
        readBinaryData(_m,input);
        readVector(_value,input);
        readVector(_rowStart,input);
		return input.good();
    }
    //mat vec
	template <typename VEC>
	void multiplyRow(const VEC& x){
		for(sizeType i=0; i<_n; ++i)
            for(sizeType j=_rowStart[i]; j<_rowStart[i+1]; ++j)
                _value[j].first*=x[i];
	}
    template <typename VEC,typename VEC_OUT>
    void multiply(const VEC& x,VEC_OUT& result) const {
        //ASSERT(_m==x.size() && _n==result.size());
        //result.resize(_n);
		#pragma omp parallel for
        for(sizeType i=0; i<_n; ++i) {
            result[i]=(T)0.0f;
            for(sizeType j=_rowStart[i]; j<_rowStart[i+1]; ++j)
                result[i]+=_value[j].first*x[_value[j].second];
        }
    }
    template <typename VEC,typename VEC_OUT>
    void multiplySubtract(const VEC& x,VEC_OUT& result) const {
        //ASSERT(_m==x.size() && _n==result.size());
        //result.resize(_n);
		#pragma omp parallel for
        for(sizeType i=0; i<_n; ++i) {
            for(sizeType j=_rowStart[i]; j<_rowStart[i+1]; ++j)
                result[i]-=_value[j].first*x[_value[j].second];
        }
    }
    template <typename VEC,typename VEC_OUT>
    void multiplyAdd(const VEC& x,VEC_OUT& result) const {
        //ASSERT(_m==x.size() && _n==result.size());
        //result.resize(_n);
		#pragma omp parallel for
        for(sizeType i=0; i<_n; ++i) {
            for(sizeType j=_rowStart[i]; j<_rowStart[i+1]; ++j)
                result[i]+=_value[j].first*x[_value[j].second];
        }
    }
    //mat vec transposed
    template <typename VEC,typename VEC_OUT>
    void multiplyTranspose(const VEC& x,VEC_OUT& result) const {
        //ASSERT(_n==x.size());
        //result.resize(_m);
        OMP_PARALLEL_FOR_
        for(sizeType i=0; i<_m; ++i)
            result[i]=(T)0.0f;
        multiplyTransposeAdd(x,result);
    }
    template <typename VEC,typename VEC_OUT>
    void multiplyTransposeSubtract(const VEC& x,VEC_OUT& result) const {
        //ASSERT(_n==x.size());
        //result.resize(_m);
        for(sizeType i=0; i<_n; ++i) {
            for(sizeType j=_rowStart[i]; j<_rowStart[i+1]; ++j)
                result[_value[j].second]-=_value[j].first*x[i];
        }
    }
    template <typename VEC,typename VEC_OUT>
    void multiplyTransposeAdd(const VEC& x,VEC_OUT& result) const {
        //ASSERT(_n==x.size());
        //result.resize(_m);
        for(sizeType i=0; i<_n; ++i) {
            for(sizeType j=_rowStart[i]; j<_rowStart[i+1]; ++j)
                result[_value[j].second]+=_value[j].first*x[i];
        }
    }
	//mat vec
    template <int B,typename VEC,typename VEC_OUT>
    void multiplyB(const VEC& x,VEC_OUT& result) const {
		sizeType i,j;
		Eigen::Matrix<Scalar,B,1> blk;
		#pragma omp parallel for private(i,j,blk)
		for(i=0; i<_n; ++i) {
			blk.setZero();
            for(j=_rowStart[i]; j<_rowStart[i+1]; ++j)
                blk+=_value[j].first*x.block<B,1>(_value[j].second*B,0);
			result.block<B,1>(i*B,0)=blk;
        }
    }
    template <int B,typename VEC,typename VEC_OUT>
    void multiplySubtractB(const VEC& x,VEC_OUT& result) const {
        for(sizeType i=0; i<_n; ++i) {
			for(sizeType j=_rowStart[i]; j<_rowStart[i+1]; ++j)
                result.block<B,1>(i*B,0)-=_value[j].first*x.block<B,1>(_value[j].second*B,0);
        }
    }
    template <int B,typename VEC,typename VEC_OUT>
    void multiplyAddB(const VEC& x,VEC_OUT& result) const {
        for(sizeType i=0; i<_n; ++i) {
			for(sizeType j=_rowStart[i]; j<_rowStart[i+1]; ++j)
                result.block<B,1>(i*B,0)+=_value[j].first*x.block<B,1>(_value[j].second*B,0);
        }
    }
    //mat vec transposed
    template <int B,typename VEC,typename VEC_OUT>
    void multiplyTransposeB(const VEC& x,VEC_OUT& result) const {
        result.setZero();
        multiplyTransposeAddB<B,VEC,VEC_OUT>(x,result);
    }
    template <int B,typename VEC,typename VEC_OUT>
    void multiplyTransposeSubtractB(const VEC& x,VEC_OUT& result) const {
        for(sizeType i=0; i<_n; ++i) {
            for(sizeType j=_rowStart[i]; j<_rowStart[i+1]; ++j)
                KERNEL_TYPE::subMulT(result.block<B,1>(_value[j].second*B,0),_value[j].first,x.block<B,1>(i*B,0));
        }
    }
    template <int B,typename VEC,typename VEC_OUT>
    void multiplyTransposeAddB(const VEC& x,VEC_OUT& result) const {
        for(sizeType i=0; i<_n; ++i) {
            for(sizeType j=_rowStart[i]; j<_rowStart[i+1]; ++j)
				KERNEL_TYPE::addMulT(result.block<B,1>(_value[j].second*B,0),_value[j].first,x.block<B,1>(i*B,0));
        }
    }
    //tranpose
    void transpose(FixedSparseMatrix& b) const {
        b.resize(_m,_n);
        b._rowStart.assign(_m+1,0);
        sizeType nnz=(sizeType)_value.size();
        for(sizeType i=0; i<nnz; i++)
            b._rowStart[_value[i].second+1]++;
        for(sizeType k=1; k<_m+1; k++)
            b._rowStart[k]+=b._rowStart[k-1];

        std::vector<sizeType> tmp=b._rowStart;
        b._value.resize(_value.size());
        for(sizeType i=0; i<_n; i++)
            for(sizeType j=_rowStart[i]; j<_rowStart[i+1]; j++) {
                sizeType& entry=tmp[_value[j].second];
                b._value[entry]=std::make_pair(_value[j].first,i);
                entry++;
            }
    }
    //mat mat transposed
	template <typename A,typename KERNEL_TYPEA,typename B,typename KERNEL_TYPEB>
    void multiplyMatrix(const FixedSparseMatrix<A,KERNEL_TYPEA>& b,FixedSparseMatrix<B,KERNEL_TYPEB>& out) const {
        //estimate nonzeros
        out.resize(_n,b.cols());
        out.getRowOffset()[0]=0;
        out.getValue().clear();
        out.getValue().reserve(_value.size()+b.getValue().size());

        //temporary rows
        sizeType nnz=0;
        sizeType rowNnz=0;
        std::vector<sizeType> off(b.cols(),-1);
        std::vector<sizeType> col(b.cols(),-1);
        for(sizeType i=0; i<_n; i++) {
            for(sizeType j=_rowStart[i]; j<_rowStart[i+1]; j++) {
                ROWE cv=_value[j];
                for(sizeType k=b.getRowOffset()[cv.second]; k<b.getRowOffset()[cv.second+1]; k++) {
					typename FixedSparseMatrix<A,KERNEL_TYPEA>::ROWE cvOut=b.getValue()[k];
                    if(off[cvOut.second] == -1) {
                        off[cvOut.second]=nnz;
                        col[rowNnz++]=cvOut.second;
                        out.getValue().push_back( std::make_pair<B,sizeType>(B(cv.first*b.getValue()[k].first), sizeType(cvOut.second)) );
                        nnz++;
                    } else {
                        out.getValue()[off[cvOut.second]].first+=cv.first*b.getValue()[k].first;
                    }
                }
            }

            out.getRowOffset()[i+1]=nnz;
            for(sizeType j=0; j<rowNnz; j++)
                off[col[j]]=-1;
            rowNnz=0;
            std::sort(out.getValue().begin()+out.getRowOffset()[i],
                      out.getValue().begin()+out.getRowOffset()[i+1],
                      typename SparseMatrix<B,KERNEL_TYPEB>::SortByCol());
        }
    }
    template <typename A,typename KERNEL_TYPEA,typename B,typename KERNEL_TYPEB>
	void multiplyMatrixTransposed(const FixedSparseMatrix<A,KERNEL_TYPEA>& b,FixedSparseMatrix<B,KERNEL_TYPEB>& out) const {
        FixedSparseMatrix<A,KERNEL_TYPEA> bTmp;
        b.transpose(bTmp);
        multiplyMatrix(bTmp,out);
    }
    template <typename A,typename KERNEL_TYPEA,typename B,typename KERNEL_TYPEB>
	void transposedMultiplyMatrix(const FixedSparseMatrix<A,KERNEL_TYPEA>& b,FixedSparseMatrix<B,KERNEL_TYPEB>& out) const {
        FixedSparseMatrix ATmp;
        transpose(ATmp);
        ATmp.multiplyMatrix(b,out);
    }
    //add sub
    void add(const FixedSparseMatrix& other,const T& coef=1.0f,const T& eps=ScalarUtil<T>::scalar_eps) {
        ASSERT(other._n == _n && other._m == _m);

        typedef typename ROW::const_iterator CITER;

        std::vector<sizeType> rowStart(_n+1,0);
        ROW value;
        value.reserve(_value.size()+other._value.size());
        std::back_insert_iterator<ROW> begV(value);
        for(sizeType i=0; i<_n; i++) {
            CITER beg=_value.begin()+_rowStart[i];
            CITER end=_value.begin()+_rowStart[i+1];
            CITER begO=other._value.begin()+other._rowStart[i];
            CITER endO=other._value.begin()+other._rowStart[i+1];
            while(beg!=end && begO!=endO) {
                if(beg->second < begO->second) {
                    *begV=*beg;
                    begV++;
                    beg++;
                } else if(beg->second > begO->second) {
                    *begV=std::make_pair(begO->first*coef,begO->second);
                    begV++;
                    begO++;
                } else {
                    T V=beg->first+begO->first*coef;
                    if(std::abs(V) > eps) {
                        *begV=std::make_pair(V,begO->second);
                        begV++;
                    }
                    begO++;
                    beg++;
                }
            }
            while(beg!=end) {
                *begV=*beg;
                begV++;
                beg++;
            }
            while(begO!=endO) {
                *begV=std::make_pair(begO->first*coef,begO->second);
                begV++;
                begO++;
            }
            rowStart[i+1]=(sizeType)value.size();
        }
        _value.swap(value);
        _rowStart.swap(rowStart);
    }
    void sub(const FixedSparseMatrix& other,const T& eps=ScalarUtil<T>::scalar_eps) {
        add(other,-1.0f,eps);
    }
	template <typename A>
	void mul(const A& coef) {
		OMP_PARALLEL_FOR_
		for(sizeType i=0;i<(sizeType)_value.size();i++)
			_value[i].first*=coef;
	}
    //indexer
    T operator()(sizeType i, sizeType j) const {
        typedef typename ROW::const_iterator ITER;
        ITER beg=_value.begin()+_rowStart[i],end=_value.begin()+_rowStart[i+1];
        ITER lb=std::lower_bound(beg,end,ROWE(KERNEL_TYPE::ZeroV(),j),typename SparseMatrix<T,KERNEL_TYPE>::SortByCol());
        if(lb == end || lb->second != j)
            return KERNEL_TYPE::ZeroV();
        else return lb->first;
    }
    void setElement(sizeType i,sizeType j,T newValue) {
        typedef typename ROW::iterator ITER;
        ROWE iv(newValue,j);
        ITER beg=_value.begin()+_rowStart[i],end=_value.begin()+_rowStart[i+1];
        ITER lb=std::lower_bound(beg,end,iv,typename SparseMatrix<T,KERNEL_TYPE>::SortByCol());
        if(lb != end && lb->second == j)
            lb->first=newValue;
        else {
            _value.insert(lb,iv);
            for(sizeType j=i+1; j<_n+1; j++)
                _rowStart[j]++;
        }
    }
    void addToElement(sizeType i,sizeType j,T incValue) {
        typedef typename ROW::iterator ITER;
        ROWE iv(incValue,j);
        ITER beg=_value.begin()+_rowStart[i],end=_value.begin()+_rowStart[i+1];
        ITER lb=std::lower_bound(beg,end,iv,typename SparseMatrix<T,KERNEL_TYPE>::SortByCol());
        if(lb != end && lb->second == j)
            lb->first+=incValue;
        else {
            _value.insert(lb,iv);
            for(sizeType j=i+1; j<_n+1; j++)
                _rowStart[j]++;
        }
    }
    //add row
    void addSparseRow(sizeType i,const ROW &r,const T& eps=ScalarUtil<T>::scalar_eps) {
        ROW newR;
        SparseMatrix<T,KERNEL_TYPE>::merge(r.begin(),r.end(),_value.begin()+_rowStart[i],_value.begin()+_rowStart[i+1],newR);

        sizeType nr=(sizeType)newR.size()-(_rowStart[i+1]-_rowStart[i]);
        if(nr > 0) {
            _value.resize(_value.size()+nr);
            for(sizeType j=_rowStart.back()-1; j>=_rowStart[i+1]; j--)
                _value[j+nr]=_value[j];
        } else if(nr < 0) {
            for(sizeType j=_rowStart[i+1]; j<_rowStart.back(); j++)
                _value[j+nr]=_value[j];
            _value.resize(_value.size()+nr);
        }
        for(sizeType j=i+1; j<_n+1; j++)
            _rowStart[j]+=nr;
        std::copy(newR.begin(),newR.end(),_value.begin()+_rowStart[i]);
    }
    void addSparseRow(sizeType i,const std::vector<sizeType>& index,const std::vector<T>& value,const T& eps=ScalarUtil<T>::scalar_eps) {
        ROW tmp;
        for(sizeType j=0; j<(sizeType)index.size(); j++)
            tmp.push_back(std::make_pair(value[j],index[j]));
        addSparseRow(i,tmp,eps);
    }
    void addSparseRowC(sizeType i,const ROW& r,const T& coef,const T& eps=ScalarUtil<T>::scalar_eps) {
        ROW valuesM=r;
        for(typename ROW::iterator beg=valuesM.begin(),end=valuesM.end(); beg!=end; beg++)
            (*beg).first*=coef;
        addSparseRow(i,valuesM,eps);
    }
    void addSparseRowC(sizeType i,const std::vector<sizeType>& index,const std::vector<T>& value,const T& coef,const T& eps=ScalarUtil<T>::scalar_eps) {
        std::vector<T> valueM=value;
        for(typename std::vector<T>::iterator beg=valueM.begin(),end=valueM.end(); beg!=end; beg++)
            (*beg)*=coef;
        addSparseRow(i,index,valueM,eps);
    }
    void addSparseRowC(sizeType i,ConstSMIterator<T> beg,ConstSMIterator<T> end,const T& coef,const T& eps=ScalarUtil<T>::scalar_eps) {
        ROW valuesM;
        for(; beg!=end; ++beg)
            valuesM.push_back(std::make_pair((*beg)*coef,beg.col()));
        addSparseRow(i,valuesM,eps);
    }
    //test
    void symmetricRemoveRolCol(sizeType i) {
        for(sizeType j=_rowStart[i]; j<_rowStart[i+1]; j++) {
            sizeType c=_value[j].second;
            _value[j].second=-1;

            typename ROW::iterator iter=
                std::lower_bound(_value.begin()+_rowStart[c],_value.begin()+_rowStart[c+1],
                                 std::pair<T,sizeType>((T)0.0f,i),SparseMatrix<T,KERNEL_TYPE>::SortByCol());
            if(iter != _value.begin()+_rowStart[c+1] && iter->second == i)
                iter->second=-1;
        }
        compact(-1);
    }
    bool isSymmetric(const T& eps=ScalarUtil<T>::scalar_eps) const {
        for(sizeType i=0; i<_n; i++)
            for(sizeType j=_rowStart[i]; j<_rowStart[i+1]; j++)
                if(std::abs(_value[j].first-operator()(_value[j].second,i)) > eps)
                    return false;
        return true;
    }
    bool isDiagonalDominant(const T& eps=ScalarUtil<T>::scalar_eps) const {
        for(sizeType i=0; i<_n; i++) {
            T res=(T)0.0f;
            for(sizeType j=_rowStart[i]; j<_rowStart[i+1]; j++)
                if(_value[j].second == i)
                    res+=std::abs(_value[j].first);
                else res-=std::abs(_value[j].first);
            if(res < -eps)
                return false;
        }
        return true;
    }
    bool equals(const FixedSparseMatrix& other,const T& eps=ScalarUtil<T>::scalar_eps) const {
        if(other.rows() != _n || other.cols() != _m)
            return false;
        for(sizeType i=0; i<_n; i++) {
            //ASSERT(numElement(i) == other.numElement(i));
            ConstSMIterator<T> cBeg=begin(i);
            ConstSMIterator<T> cEnd=end(i);
            ConstSMIterator<T> oBeg=other.begin(i);
            ConstSMIterator<T> oEnd=other.end(i);
            if(!SparseMatrix<T,KERNEL_TYPE>::rowEquals(cBeg,cEnd,oBeg,oEnd,eps))
                return false;
        }
        return true;
    }
    bool equals(const SparseMatrix<T,KERNEL_TYPE>& other,const T& eps=ScalarUtil<T>::scalar_eps) const {
        if(other.rows() != _n || other.cols() != _m)
            return false;
        for(sizeType i=0; i<_n; i++) {
            //ASSERT(numElement(i) == other.numElement(i));
            ConstSMIterator<T> cBeg=begin(i);
            ConstSMIterator<T> cEnd=end(i);
            ConstSMIterator<T> oBeg=other.begin(i);
            ConstSMIterator<T> oEnd=other.end(i);
            if(!SparseMatrix<T,KERNEL_TYPE>::rowEquals(cBeg,cEnd,oBeg,oEnd,eps))
                return false;
        }
        return true;
    }
    bool isZero(const T& eps=ScalarUtil<T>::scalar_eps) const {
        sizeType nnz=(sizeType)_value.size();
        for(sizeType i=0; i<nnz; i++)
            if(std::abs(_value[i].first) > eps)
                return false;
        return true;
    }
    void compact(sizeType flag) {
        sizeType i=0;
        sizeType j=0;
        sizeType k=0;
        sizeType nnz=(sizeType)_value.size();
        for(; i<nnz; i++) {
            while(i == _rowStart[k+1] && k<_n) {
                _rowStart[k+1]=j;
                k++;
            }
            if(_value[i].second != flag) {
                _value[j]=_value[i];
                j++;
            }
        }
        _value.resize(j);
        for(; k<_n; k++)
            _rowStart[k+1]=j;
    }
    //iterator
    ConstSMIterator<T> begin(sizeType r) const {
        return ConstSMIterator<T>(r,_value.begin()+_rowStart[r]);
    }
    ConstSMIterator<T> end(sizeType r) const {
        return ConstSMIterator<T>(r,_value.begin()+_rowStart[r+1]);
    }
    SMIterator<T> begin(sizeType r) {
        return SMIterator<T>(r,_value.begin()+_rowStart[r]);
    }
    SMIterator<T> end(sizeType r) {
        return SMIterator<T>(r,_value.begin()+_rowStart[r+1]);
    }
    //dangerous method for optimized routines only
    FORCE_INLINE const ROW& getValue() const{return _value;}
    FORCE_INLINE ROW& getValue(){return _value;}
	FORCE_INLINE const std::vector<sizeType>& getRowOffset() const{return _rowStart;}
	FORCE_INLINE std::vector<sizeType>& getRowOffset(){return _rowStart;}
    void getMemoryConsumption(T& MBA) const
    {
        MBA=(T)(sizeof(ROWE)*_value.size()+sizeof(sizeType)*_rowStart.size());
        MBA/=(1024.0f*1024.0f);
    }
	void printFillin() const
	{
		INFOV("Fillin: %f",(scalarD)_value.size()/(scalarD)rows())
	}
protected:
    sizeType _n;
    sizeType _m;
    ROW _value;
    std::vector<sizeType> _rowStart;
};

PRJ_END

#endif
