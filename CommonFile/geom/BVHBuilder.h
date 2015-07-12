#ifndef BVH_BUILDER_H
#define BVH_BUILDER_H

#include "MathBasic.h"
#include "IO.h"
#include <deque>
#include <stack>
#include <boost/static_assert.hpp>

USE_PRJ_NAMESPACE

namespace GEOM{

  template <int D>struct SurfaceArea;
  template <>struct SurfaceArea<3> {
    static scalar area(const BBox<scalar>& bb) {
	  Vec3 ext=compMax(Vec3::Zero(),(Vec3)(bb._maxC-bb._minC));
	  return (ext[0]*ext[1]+ext[0]*ext[2]+ext[1]*ext[2])*2.0f;
    }
  };
  template <>struct SurfaceArea<2> {
    static scalar area(const BBox<scalar>& bb) {
	  Vec3 ext=compMax(Vec3::Zero(),(Vec3)(bb._maxC-bb._minC));
	  return (ext[0]+ext[1])*2.0f;
    }
  };
  template <typename NODE,int DIM>
  struct BVHBuilder {
  public:
    BOOST_STATIC_ASSERT(DIM == 2 || DIM == 3);
    struct BVHHandle {
	  BVHHandle():_left(false),_cost(0.0f) {}
	  BVHHandle(scalar val,sizeType rid,bool left):_val(val),_cost(0.0f),_rid(rid),_left(left) {}
	  bool operator<(const BVHHandle& other) const {
		return _val<other._val;
	  }
	  bool operator<=(const BVHHandle& other) const {
		return _val<=other._val;
	  }
	  bool operator>(const BVHHandle& other) const {
		return _val>other._val;
	  }
	  bool operator>=(const BVHHandle& other) const {
		return _val>=other._val;
	  }
	  bool operator==(const BVHHandle& other) const {
		return _val==other._val;
	  }
	  scalar _val,_cost;
	  sizeType _rid;
	  bool _left;
    };
    sizeType buildBVH(vector<NODE>& bvh) {
	  if(bvh.empty())
		return -1;
	  //find roots
	  _roots.clear();
	  for(sizeType i=0; i<(sizeType)bvh.size(); i++)
		if(bvh[i]._parent == -1)
		  _roots.push_back(i);
	  //short bounding box alone each axis
	  for(sizeType d=0; d<DIM; d++) {
		for(sizeType i=0; i<(sizeType)_roots.size(); i++) {
		  _hdls[d].push_back(BVHHandle(bvh[_roots[i]]._bb.minCorner()[d],_roots[i],true));
		  _hdls[d].push_back(BVHHandle(bvh[_roots[i]]._bb.maxCorner()[d],_roots[i],false));
		}
		std::sort(_hdls[d].begin(),_hdls[d].end());
	  }
	  //build BVH
	  Vec3i f(0,0,0);
	  Vec3i t((sizeType)_hdls[0].size(),
			  (sizeType)_hdls[1].size(),
			  (sizeType)_hdls[2].size());
	  bvh.reserve(bvh.size()+_roots.size()-1);
	  sizeType ret=buildBVH(bvh,0,(sizeType)_roots.size(),f,t);
	  bvh.back()._parent=-1;
	  return ret;
    }
  private:
	//#define DEBUG_BVH
    sizeType buildBVH(vector<NODE>& bvh,sizeType fr,sizeType tr,Vec3i f,Vec3i t) {
	  if(tr == fr+1)
		return _roots[fr];
	  //find split location
	  sizeType lt=fr,rf=tr;
	  Vec3i ltd=f,rfd=t;
	  {
		//calculate split cost
		sizeType minD=-1,minId=-1;
		scalar minCost=numeric_limits<scalar>::max();
		for(sizeType d=0; d<DIM; d++) {
#ifdef DEBUG_BVH
		  {
			std::vector<BVHHandle> hdls;
			for(sizeType i=fr; i<tr; i++) {
			  hdls.push_back(BVHHandle(bvh[_roots[i]]._bb._minC[d],_roots[i],true));
			  hdls.push_back(BVHHandle(bvh[_roots[i]]._bb._maxC[d],_roots[i],false));
			}
			std::sort(hdls.begin(),hdls.end());

			ASSERT((sizeType)hdls.size() == t[d]-f[d])
			  ASSERT((sizeType)hdls.size() == (tr-fr)*2)
			  for(sizeType i=0; i<(sizeType)hdls.size(); i++)
				ASSERT(hdls[i] == _hdls[d][i+f[d]]);
		  }
#endif
		  calcCost(bvh,_hdls[d],f[d],t[d]);
		  for(sizeType i=f[d]; i<(sizeType)t[d]; i++) {
#ifdef DEBUG_BVH
			ASSERT(debugCost(bvh,_hdls[d],f[d],t[d],i) == _hdls[d][i]._cost);
#endif
			if(_hdls[d][i]._cost < minCost) {
			  minCost=_hdls[d][i]._cost;
			  minD=d;
			  minId=i;
			}
		  }
		}
		//find split set
		for(sizeType i=fr; i<tr; i++)
		  bvh[_roots[i]]._parent=-1;
		const BVHHandle& minHdl=_hdls[minD][minId];
		for(sizeType i=f[minD]; i<t[minD]; i++) {
		  const BVHHandle& hdl=_hdls[minD][i];
		  if(bvh[hdl._rid]._parent != -1)
			continue;
		  if(hdl._left && hdl <= minHdl) {
			_roots[lt++]=hdl._rid;
			bvh[hdl._rid]._parent=1;
		  }
		  if(!hdl._left && hdl >= minHdl) {
			_roots[--rf]=hdl._rid;
			bvh[hdl._rid]._parent=2;
		  }
		}
		ASSERT_MSG(lt == rf,"Strange Error!");
		if(lt == tr || rf == fr) {
		  lt=rf=(fr+tr)/2;
		  for(sizeType i=fr; i<lt; i++)
			bvh[_roots[i]]._parent=1;
		  for(sizeType i=rf; i<tr; i++)
			bvh[_roots[i]]._parent=2;
		}
		for(sizeType d=0; d<DIM; d++)
		  split(bvh,_hdls[d],ltd[d],rfd[d]);
	  }

	  //split
	  NODE n;
	  n._l=buildBVH(bvh,fr,lt,f,ltd);
	  n._r=buildBVH(bvh,rf,tr,rfd,t);
	  bvh[n._l]._parent=(sizeType)bvh.size();
	  bvh[n._r]._parent=(sizeType)bvh.size();

	  //add new root
	  n._bb=bvh[n._l]._bb;
	  n._bb.setUnion(bvh[n._r]._bb);
	  n._nrCell=bvh[n._l]._nrCell+bvh[n._r]._nrCell;
	  bvh.push_back(n);
	  return bvh.size()-1;
    }
    static void calcCost(const vector<NODE>& bvh,vector<BVHHandle>& X,sizeType f,sizeType t) {
	  BBox<scalar> bb;
	  sizeType nrC=0;
	  for(sizeType i=f; i<t; i++) {
		X[i]._cost=SurfaceArea<DIM>::area(bb)*(scalar)nrC;
		if(X[i]._left) {
		  Vec3 minC=bvh[X[i]._rid]._bb.minCorner();
		  Vec3 maxC=bvh[X[i]._rid]._bb.maxCorner();
		  bb.setUnion(BBox<scalar>(minC,maxC));
		  nrC+=bvh[X[i]._rid]._nrCell;
		}
	  }
	  bb.reset();
	  nrC=0;
	  for(sizeType i=t-1; i>=f; i--) {
		X[i]._cost+=SurfaceArea<DIM>::area(bb)*(scalar)nrC;
		if(!X[i]._left) {
		  Vec3 minC=bvh[X[i]._rid]._bb.minCorner();
		  Vec3 maxC=bvh[X[i]._rid]._bb.maxCorner();
		  bb.setUnion(BBox<scalar>(minC,maxC));
		  nrC+=bvh[X[i]._rid]._nrCell;
		}
	  }
    }
    static scalar debugCost(const vector<NODE>& bvh,vector<BVHHandle>& X,sizeType f,sizeType t,sizeType I) {
	  scalar cost=0.0f;
	  BBox<scalar> bb;
	  sizeType nrC=0;
	  for(sizeType i=f; i<I; i++)
		if(X[i]._left) {
		  bb.setUnion(bvh[X[i]._rid]._bb);
		  nrC+=bvh[X[i]._rid]._nrCell;
		}
	  cost+=SurfaceArea<DIM>::area(bb)*(scalar)nrC;

	  bb.reset();
	  nrC=0;
	  for(sizeType i=t-1; i>I; i--)
		if(!X[i]._left) {
		  bb.setUnion(bvh[X[i]._rid]._bb);
		  nrC+=bvh[X[i]._rid]._nrCell;
		}
	  cost+=SurfaceArea<DIM>::area(bb)*(scalar)nrC;

	  return cost;
    }
    static void split(const vector<NODE>& bvh,vector<BVHHandle>& X,sizeType& lt,sizeType& rf) {
	  std::vector<BVHHandle> R;
	  sizeType j=lt;
	  for(sizeType i=lt; i<rf; i++) {
		sizeType p=bvh[X[i]._rid]._parent;
		ASSERT(p == 1 || p == 2)
		  if(p == 1) X[j++]=X[i];
		  else R.push_back(X[i]);
	  }
	  lt=rf=j;
	  std::copy(R.begin(),R.end(),X.begin()+rf);
    }
    std::vector<BVHHandle> _hdls[3];
    std::vector<sizeType> _roots;
  };

  template <typename T,typename BBOX=BBox<scalar> >
  struct Node : public Serializable {
    Node():Serializable(4),_l(-1),_r(-1),_parent(-1),_nrCell(-1) {}
    bool read(std::istream& is,IOData* dat) {
	  readBinaryData(_bb,is);
	  readBinaryData(_cell,is,dat);
	  readBinaryData(_l,is);
	  readBinaryData(_r,is);
	  readBinaryData(_parent,is);
	  readBinaryData(_nrCell,is);
	  return is.good();
    }
    bool write(std::ostream& os,IOData* dat) const {
	  writeBinaryData(_bb,os);
	  writeBinaryData(_cell,os,dat);
	  writeBinaryData(_l,os);
	  writeBinaryData(_r,os);
	  writeBinaryData(_parent,os);
	  writeBinaryData(_nrCell,os);
	  return os.good();
    }
    boost::shared_ptr<Serializable> copy() const {
	  return boost::shared_ptr<Serializable>(new Node<T,BBOX>);
    }
    BBOX _bb;
    T _cell;
    sizeType _l,_r,_parent,_nrCell;
  };
  template <typename T,typename BBOX>
  void buildBVH(vector<Node<T,BBOX> >& bvh,sizeType dim,T verbose)
  {
    sizeType nrBVH=(sizeType)bvh.size();
    if(dim == 2)BVHBuilder<Node<T,BBOX>,2>().buildBVH(bvh);
    else BVHBuilder<Node<T,BBOX>,3>().buildBVH(bvh);
    for(sizeType i=nrBVH; i<(sizeType)bvh.size(); i++)
	  bvh[i]._cell=verbose;
  }

  template <typename T,typename BBOX>
  static void writeBVHByLevel(const vector<Node<T,BBOX> >& bvh,T verbose)
  {
    if(bvh.empty())
	  return;
    deque<sizeType> lv;
    lv.push_back((sizeType)bvh.size()-1);
    for(sizeType i=0; !lv.empty(); i++) {
	  ostringstream oss;
	  oss << "./lv" << i << ".vtk";
	  VTKWriter<scalar> os("BB Level",oss.str(),true);

	  sizeType nrN=(sizeType)lv.size();
	  std::vector<Vec3,Eigen::aligned_allocator<Vec3> > bbs;
	  for(sizeType i=0; i<nrN; i++) {
		const Node<T,BBOX>& n=bvh[lv.front()];
		BBOX bb=n._bb;
		bbs.push_back(bb.minCorner());
		bbs.push_back(bb.maxCorner());
		if(n._cell == verbose) {
		  lv.push_back(n._l);
		  lv.push_back(n._r);
		}
		lv.pop_front();
	  }
	  os.appendVoxels(bbs.begin(),bbs.end(),false);
    }
  }
  template <typename T,typename BBOX=BBox<scalar> >
  class BVHQuery
  {
  public:
    template <typename T2,typename CACHE>
    struct AddCache {
	  AddCache(std::vector<CACHE>& cache):_cache(cache) {}
	  void onCell(const Node<T,BBOX>& n,const Node<T2,BBOX>& n2) {
		_cache.push_back(CACHE(n._cell,n2._cell));
	  }
	  std::vector<CACHE>& _cache;
    };
    BVHQuery(const vector<Node<T,BBOX> >& bvh,sizeType dim,T verbose)
	  :_bvh(bvh),_dim(dim),_verbose(verbose),_active(NULL) {}
    void updateBVH(scalar expand=1E-3f) {
	  sizeType nrN=(sizeType)_bvh.size();
	  for(sizeType i=0; i<nrN; i++) {
		Node<T,BBOX>& n=const_cast<Node<T,BBOX>&>(_bvh[i]);
		if(n._cell != _verbose)
		  n._bb.enlarged(expand,_dim);
		else {
		  n._bb=_bvh[n._l]._bb;
		  n._bb.setUnion(_bvh[n._r]._bb);
		}
	  }
    }
    template <typename CALLBACK>
    void pointQuery(const Vec3& pos,scalar eps,CALLBACK& cb,scalar& sqrDist) const {
	  if(_bvh.empty())
		return;
	  std::stack<sizeType> ss;
	  ss.push((sizeType)_bvh.size()-1);
	  while(!ss.empty()) {
		const Node<T,BBOX>& node=_bvh[ss.top()];
		ss.pop();
		if(node._bb.distTo(pos,_dim) > eps)
		  continue;
		if(node._cell == _verbose) {
		  ss.push(node._l);
		  ss.push(node._r);
		} else {
		  cb.updateDist(node,sqrDist);
		  if(sqrDist == 0.0f)return;
		}
	  }
    }
    template <typename CALLBACK>
    void pointQuery(CALLBACK& cb) const {
	  if(_bvh.empty())
		return;
	  std::stack<sizeType> ss;
	  ss.push((sizeType)_bvh.size()-1);
	  while(!ss.empty()) {
		const Node<T,BBOX>& node=_bvh[ss.top()];
		ss.pop();
		if(!cb.validNode(node))
		  continue;
		if(node._cell == _verbose) {
		  ss.push(node._l);
		  ss.push(node._r);
		} else {
		  cb.updateDist(node);
		}
	  }
    }
    template <typename CALLBACK>
    void pointDistQuery(const Vec3& pt,CALLBACK& cb,Vec3& cp,Vec3& n,scalar& dist) const {
	  if(_bvh.empty())
		return;
	  std::stack<sizeType> ss;
	  ss.push((sizeType)_bvh.size()-1);
	  while(!ss.empty()) {
		const Node<T,BBOX>& node=_bvh[ss.top()];
		ss.pop();
		if(node._bb.distTo(pt,_dim) > std::min<scalar>(cb.depth(),dist))
		  continue;
		if(node._cell == _verbose) {
		  ss.push(node._l);
		  ss.push(node._r);
		} else {
		  cb.updateDist(node,pt,cp,n,dist);
		}
	  }
    }
    template <typename T2,typename CALLBACK>
    void interBodyQuery(const BVHQuery<T2,BBOX>& query2,CALLBACK& cb) const {
	  if(_bvh.empty() || query2._bvh.empty())
		return;
	  std::stack<Vec2i> ss;
	  ss.push(Vec2i(_bvh.size()-1,query2._bvh.size()-1));
	  while(!ss.empty()) {
		sizeType rA=ss.top()[0];
		sizeType rB=ss.top()[1];
		ss.pop();

		if(_active && query2._active && !(*_active)[rA] && !(*(query2._active))[rB])continue;
		if(!_bvh[rA]._bb.intersect(query2._bvh[rB]._bb,_dim))continue;
		if(_bvh[rA]._cell != _verbose && query2._bvh[rB]._cell != query2._verbose) {
		  cb.onCell(_bvh[rA],query2._bvh[rB]);
		} else if(_bvh[rA]._cell != _verbose && query2._bvh[rB]._cell == query2._verbose) {
		  ss.push(Vec2i(rA,query2._bvh[rB]._l));
		  ss.push(Vec2i(rA,query2._bvh[rB]._r));
		} else if(_bvh[rA]._cell == _verbose && query2._bvh[rB]._cell != query2._verbose) {
		  ss.push(Vec2i(_bvh[rA]._l,rB));
		  ss.push(Vec2i(_bvh[rA]._r,rB));
		} else {
		  ss.push(Vec2i(_bvh[rA]._l,query2._bvh[rB]._l));
		  ss.push(Vec2i(_bvh[rA]._l,query2._bvh[rB]._r));
		  ss.push(Vec2i(_bvh[rA]._r,query2._bvh[rB]._l));
		  ss.push(Vec2i(_bvh[rA]._r,query2._bvh[rB]._r));
		}
	  }
    }
    template <typename T2,typename CALLBACK>
    void interBodyDepthQuery(const BVHQuery<T2,BBOX>& query2,CALLBACK& cb,scalar depth) const {
	  if(_bvh.empty() || query2._bvh.empty())
		return;
	  std::stack<Vec2i> ss;
	  ss.push(Vec2i(_bvh.size()-1,query2._bvh.size()-1));
	  while(!ss.empty()) {
		sizeType rA=ss.top()[0];
		sizeType rB=ss.top()[1];
		ss.pop();

		if(_bvh[rA]._bb.distTo(query2._bvh[rB]._bb,_dim) > depth)continue;
		if(_bvh[rA]._cell != _verbose && query2._bvh[rB]._cell != query2._verbose) {
		  cb.onCell(_bvh[rA],query2._bvh[rB]);
		} else if(_bvh[rA]._cell != _verbose && query2._bvh[rB]._cell == query2._verbose) {
		  ss.push(Vec2i(rA,query2._bvh[rB]._l));
		  ss.push(Vec2i(rA,query2._bvh[rB]._r));
		} else if(_bvh[rA]._cell == _verbose && query2._bvh[rB]._cell != query2._verbose) {
		  ss.push(Vec2i(_bvh[rA]._l,rB));
		  ss.push(Vec2i(_bvh[rA]._r,rB));
		} else {
		  ss.push(Vec2i(_bvh[rA]._l,query2._bvh[rB]._l));
		  ss.push(Vec2i(_bvh[rA]._l,query2._bvh[rB]._r));
		  ss.push(Vec2i(_bvh[rA]._r,query2._bvh[rB]._l));
		  ss.push(Vec2i(_bvh[rA]._r,query2._bvh[rB]._r));
		}
	  }
    }
    template <typename T2,typename CACHE>
    void broadphaseQuery(const BVHQuery<T2,BBOX>& query2,std::vector<CACHE>& cache) const {
	  AddCache<T2,CACHE> cb(cache);
	  interBodyQuery(query2,cb);
    }
    const vector<bool>* _active;
	const vector<Node<T,BBOX> >& _bvh;
    sizeType _dim;
    T _verbose;
  };

}

#endif
