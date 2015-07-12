#ifndef _SELECTIONGROUP_H_
#define _SELECTIONGROUP_H_

#include <vector>
#include <set>
#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
using std::vector;
using std::set;
using std::string;

namespace UTILITY{
  
  /**
   * @class SelectionGroup manager selection group, including adding, remomve
   * groups.
   */
  class SelectionGroup{
	
  public:
	void setGroup(const vector<set<int> > &groups);
	void setGroup(const vector<vector<int> > &groups);
	set<int> addGroup( const set<int>& );
	set<int> addGroup( const vector<int>& );

	set<int> removeGroup( const set<int>& );
	set<int> removeGroup( const vector<int>& );

	int removeGroup( set<int>& group)const;
	int removeGroup( vector<int>& group)const;

	bool contain(const int ele_id)const;
	bool isEmpty()const{
	  return groups.size() <= 0;
	}
	const vector<set<int> > &getGroup()const{
	  return groups;
	}
	vector<vector<int> > getGroupVec()const;
	template<class VECTOR_INT>
	void getVector(VECTOR_INT &vec)const{
	  vec.clear();
	  vec.reserve(numNodes());
	  BOOST_FOREACH(const set<int>& s, groups){
		BOOST_FOREACH(const int i, s){
		  vec.push_back(i);
		};
	  }
	}
	int numGroup()const{
	  return (int)groups.size();
	}
	int numNodes()const{
	  int count = 0;
	  for (size_t i = 0; i < groups.size(); ++i)
		count += groups[i].size();
	  return count;
	}

	void clear(){
	  groups.clear();
	}
	bool saveAsVec(const string file_name)const;
	bool load(const string file_name);

  protected:
	set<int> vec2set(const vector<int> &)const;
	vector<int> set2vec(const set<int> &)const;
	
  private:
	vector<set<int> > groups;
	
  };
  
  typedef boost::shared_ptr<SelectionGroup> pSelectionGroup;
  
}//end of namespace

#endif /*_SELECTIONGROUP_H_*/
