#ifndef DISJOINT_SET_H
#define DISJOINT_SET_H

#include "MathBasic.h"
#include "Zero.h"

PRJ_BEGIN

template <typename WEIGHT>
struct DisjointSetElem {
    sizeType _rank;	//rank based merging
    sizeType _p;	//pointer to parent
    sizeType _size;	//size of current subtree, size=0 is a tag to avoid joining
    WEIGHT _Int;	//another criterion for merging
};

template <typename WEIGHT,typename ELEM=DisjointSetElem<WEIGHT> >
class DisjointSet
{
public:
    DisjointSet():_nrSet(0) {}
    DisjointSet(sizeType elements) {
        resize(elements);
    }
    void resize(sizeType elements) {
        _elts.clear();
        _elts.resize(elements);
        _nrSet=elements;
        for(sizeType i=0; i < elements; i++) {
            _elts[i]._rank=0;
            _elts[i]._size=1;
            _elts[i]._p=i;
            _elts[i]._Int=Zero<WEIGHT>::value();
        }
    }
    virtual ~DisjointSet() {}
    sizeType find(sizeType x) {
        sizeType y=x;
        while(y != _elts[y]._p)
            y=_elts[y]._p;
        _elts[x]._p=y;
        return y;
    }
    void join(sizeType x,sizeType y) {
        if(x == y)
            return;
        if(_elts[x]._rank > _elts[y]._rank) {
            _elts[y]._p=x;
            _elts[x]._size+=_elts[y]._size;
        } else {
            _elts[x]._p=y;
            _elts[y]._size+=_elts[x]._size;
            if(_elts[x]._rank == _elts[y]._rank)
                _elts[y]._rank++;
        }
        _nrSet--;
    }
    void joinSafe(sizeType x,sizeType y) {
        if(_elts[x]._size == 0 || _elts[y]._size == 0)
            return;
        join(find(x),find(y));
    }
    sizeType size(sizeType x) const {
        return _elts[x]._size;
    }
    sizeType numSets() const {
        return _nrSet;
    }
	std::vector<ELEM,boost::fast_pool_allocator<ELEM> > _elts;
    sizeType _nrSet;
};

PRJ_END

#endif
