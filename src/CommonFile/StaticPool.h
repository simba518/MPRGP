#ifndef STATIC_POOL_H
#define STATIC_POOL_H

#include "Config.h"
#include "IO.h"
#include <vector>

PRJ_BEGIN

template <typename T,int NR>
class StaticStack{
public:
    StaticStack():_ptr(0){}
    bool empty() const{return _ptr==0;}
    sizeType size() const{return _ptr;}
    bool push(T p){
        if(_ptr == NR)
            return false;
        else{
            _pool[_ptr++]=p;
            return true;
        }
    }
    T pop(){
        ASSERT(!empty());
        return _pool[--_ptr];
    }
protected:
    T _pool[NR+1];
    sizeType _ptr;
};
template <typename T>
class StaticPool
{
public:
    struct ConstIterator {
        typedef typename T value_type;
        ConstIterator(sizeType id,const StaticPool& pool):_id(id),_pool(pool) {}
        void operator++() {_id=_pool(_id)._next;}
        bool operator==(const ConstIterator& other) const {return _id == other._id;}
        bool operator!=(const ConstIterator& other) const {return !operator==(other);}
        const value_type& operator*() const {return _pool(_id);}
        sizeType _id;
        const StaticPool& _pool;
    };
    struct Iterator : public ConstIterator {
        Iterator(sizeType id,const StaticPool& pool):ConstIterator(id,pool) {}
        value_type& operator*() {return (value_type&)(_pool(_id));}
    };
public:
    StaticPool():_nr(0),_free(-1),_head(-1) {
        _pool.resize(16);
        clear();
    }
    void write(std::ostream& os) const
    {
        writeBinaryData(_nr,os);
        writeBinaryData(_free,os);
        writeBinaryData(_head,os);
        writeVector(_pool,os);
    }
    void read(std::istream& is)
    {
        readBinaryData(_nr,is);
        readBinaryData(_free,is);
        readBinaryData(_head,is);
        readVector(_pool,is);
    }
	bool writeNonAtomic(std::ostream& os) const
	{
		writeBinaryData(nr(),os);
		for(sizeType i=usedHead(); i!=-1;) {
			writeBinaryData(i,os);
			if(!(operator()(i)).write(os))
				return false;
			i=operator()(i)._next;
		}
		return os.good();
	}
	bool readNonAtomic(std::istream& is)
	{
		clear();
		sizeType nrE;
		readBinaryData(nrE,is);
		for(sizeType i=0;i<nrE;i++) {
			sizeType id;
			readBinaryData(id,is);
			ASSERT(allocate(id));

			if(!(operator()(id)).read(is))
				return false;
		}
		return is.good();
	}
    void clear() {
        for(sizeType i=0; i<(sizeType)_pool.size(); i++) {
            _pool[i]._next=i+1;
            _pool[i]._prev=i-1;
        }
        _pool.front()._prev=-1;
        _pool.back()._next=-1;
        _free=0;
        _head=-1;
        _nr=0;
    }
    bool allocate(sizeType i){
        while((sizeType)_pool.size() <= i)
            expand();

        {
            sizeType curr=i;
            while(_pool[curr]._prev>=0)
                curr=_pool[curr]._prev;
            ASSERT(curr == _free || curr == _head)
            if(curr == _head)
                return false;
        }

        //erase
        if(_free == i)
            _free=_pool[i]._next;
        sizeType n=_pool[i]._next;
        sizeType p=_pool[i]._prev;
        if(n>=0)_pool[n]._prev=p;
        if(p>=0)_pool[p]._next=n;
        _nr++;

        _pool[i]._next=_head;
        if(_head>=0)_pool[_head]._prev=i;
        _head=i;
        return true;
    }
    sizeType allocate() {
        if(_free < 0) 
            expand();

        sizeType ret=_free;
        _free=_pool[_free]._next;
        if(_free>=0)_pool[_free]._prev=-1;
        _nr++;

        _pool[ret]._next=_head;
        if(_head>=0)_pool[_head]._prev=ret;
        _head=ret;
        return _head;
    }
    ConstIterator destroy(sizeType i) {
        ConstIterator ret(_pool[i]._next,*this);
        if(_head == i)
            _head=_pool[i]._next;

        //erase
        sizeType n=_pool[i]._next;
        sizeType p=_pool[i]._prev;
        if(n>=0)_pool[n]._prev=p;
        if(p>=0)_pool[p]._next=n;

        //insert
        _pool[i]._prev=-1;
        _pool[i]._next=_free;
        if(_free>=0)_pool[_free]._prev=i;
        _free=i;
        _nr--;
        return ret;
    }
    T& operator()(sizeType i) {
        return _pool[i];
    }
    const T& operator()(sizeType i) const {
        return _pool[i];
    }
    void parityCheck() const {
        std::vector<unsigned char> sz(_pool.size(),0);
        for(sizeType i=_free; i!=-1;) {
            ASSERT(sz[i] == 0)
            sz[i]=1;
            i=_pool[i]._next;
        }
        sizeType nr=0;
        for(sizeType i=_head; i!=-1;) {
            ASSERT(sz[i] == 0)
            sz[i]=2;
            i=_pool[i]._next;
            nr++;
        }
        for(sizeType i=0; i<(sizeType)_pool.size(); i++) {
            sizeType n=_pool[i]._next;
            sizeType p=_pool[i]._prev;
            if(n >= 0) {
                ASSERT(_pool[n]._prev==i);
                ASSERT(sz[i] == sz[n]);
            }
            if(p >= 0) {
                ASSERT(_pool[p]._next==i);
                ASSERT(sz[i] == sz[p]);
            }
        }
        if(_free >= 0)
            ASSERT(_pool[_free]._prev==-1)
            if(_head >= 0)
                ASSERT(_pool[_head]._prev==-1)
                ASSERT(nr == _nr);
        for(sizeType i=0; i<(sizeType)sz.size(); i++)
            ASSERT(sz[i] != 0)
        }
    sizeType nr() const{return _nr;}
    sizeType freeHead() const{return _free;}
    sizeType usedHead() const{return _head;}
    Iterator begin(){return Iterator(_head,*this);}
    Iterator end(){return Iterator(-1,*this);}
    ConstIterator begin() const{return ConstIterator(_head,*this);}
    ConstIterator end() const{return ConstIterator(-1,*this);}
protected:
    void expand(){
        //reallocate and copy
        std::vector<T> tmp(_pool.size()*2);
        std::copy(_pool.begin(),_pool.end(),tmp.begin());

        //insert empty
        for(sizeType i=(sizeType)_pool.size(); i<(sizeType)tmp.size(); i++) {
            tmp[i]._next=i+1;
            tmp[i]._prev=i-1;
        }
        tmp.back()._next=_free;
        tmp[_pool.size()]._prev=-1;

        if(_free>=0)
            tmp[_free]._prev=(sizeType)tmp.size()-1;
        _free=_pool.size();
        _pool.swap(tmp);
    }
    sizeType _nr;
    sizeType _free;
    sizeType _head;
    std::vector<T> _pool;
};

PRJ_END

#endif