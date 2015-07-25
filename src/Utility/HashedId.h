#ifndef _HASHEDID_H_
#define _HASHEDID_H_

#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <eigen3/Eigen/Dense>

namespace UTILITY{
  
  struct HashedId{
  public:
	HashedId();
	HashedId(int a,int b,const int& no);
	HashedId(int a,int b,int c,const int& no);
	HashedId(const Eigen::Vector3i& pos,const int& edgeLen,const int& no);
	bool operator==(const HashedId& other) const;
	int _id[4];
	int _no;
  };

  struct HashedIdHash: public boost::hash<HashedId>{

	std::size_t operator()(const HashedId& i) const{
	  const boost::hash<int> h;
	  return h(i._id[0])+h(i._id[1])+h(i._id[2])+h(i._id[3]);
	}
  };
}

#endif /* _HASHEDID_H_ */

