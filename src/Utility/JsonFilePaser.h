#ifndef _JSONFILEPASER_H_
#define _JSONFILEPASER_H_

#include <stdlib.h>
#include <vector>
#include <set>
#include <fstream>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/filesystem.hpp>
#include <Log.h>
#include <IOTools.h>
using namespace boost::property_tree;

namespace UTILITY{
  
  /**
   * @class JsonFilePaser base interface for reading the initfile.
   * 
   */
  class JsonFilePaser{
	
  public:
	JsonFilePaser(bool print = false);

	~JsonFilePaser();

	bool open(const std::string filename);
	void close();

	std::string getFileDir()const;
	std::string getFileName()const;

	template<class T>
	bool read(const std::string eleName, T &value){
	  if( !_file.is_open() ){
		ERROR_LOG("json file is not open");
		return false;
	  }
	  bool succ = true;
	  try{
		actualRead(eleName,value);
	  }catch(std::exception& ex){
		succ = false;
		ERROR_LOG(ex.what());
	  }
	  return succ;
	}
	template<class T>
	bool read(const std::string eleName, T &value,const T default_value){
	  const bool succ = read(eleName, value);
	  if(!succ){
		value = default_value;
		WARN_LOG("use default value for "<<eleName<<": "<<default_value);
	  }
	  return succ;
	}
	template<class VECTOR3D>
	bool read3d(const std::string eleName, VECTOR3D &vec){
	  assert_eq(vec.size(), 3);
	  std::vector<double> tv;
	  bool succ = false;
	  if(read(eleName, tv) && tv.size() == 3){
		vec[0] = tv[0];
		vec[1] = tv[1];
		vec[2] = tv[2];
		succ = true;
	  }
	  ERROR_LOG_COND("the length of the vector should be 3, not "<< tv.size(), tv.size()==3);
	  return succ;
	}

	template<class VECTOR>
	bool readVecFile(const std::string eleName, VECTOR &value,IO_TYPE io_type=BINARY){
	  std::string vectorFile;
	  return readFilePath(eleName,vectorFile,true)?loadVec(vectorFile,value,io_type):false;
	}
	template<class T>
	bool readVecFile(const std::string eleName, std::set<T> &value,IO_TYPE io_type=BINARY){
	  std::vector<T> v;
	  const bool succ = readVecFile(eleName, v, io_type);
	  if (succ){
		value.clear();
		for (size_t i = 0; i < v.size(); ++i)
		  value.insert(v[i]);
	  }
	  return succ;
	}
	template<class MATRIX>
	bool readMatFile(const std::string eleName, MATRIX &value,IO_TYPE io_type=BINARY){
	  std::string matFile;
	  return readFilePath(eleName,matFile,true)?loadMat(matFile,value,io_type):false;
	}
	bool readFilePath(const std::string eleName, std::string &filePath,const bool checkFileExist = true);
	bool readFilePath(const std::string eleName, std::string &filePath,const std::string default_value,const bool checkFileExist = true){
	  const bool succ = readFilePath(eleName, filePath, checkFileExist);
	  if (!succ){
		filePath = default_value;
		WARN_LOG("use default value for "<<eleName<<": "<<default_value);
	  }
	  return succ;
	}
	bool readFilePath(const std::string eleName, std::vector<std::string> &filePathes,const bool checkFileExist = true);
	std::string replaceHomeDir(const std::string &path,const bool checkPathExist = true);
	
  protected:
	bool fileIsExist(const std::string eleName,const std::string filePath,bool check=true)const;
	bool dirIsExist(const std::string eleName,const std::string dirPath,bool check)const;

	template<class T>
	void actualRead(const std::string eleName, T &value){
	  value = _jsonData.get<T>(eleName);
	}

	template<class T>
	void actualRead(const std::string eleName, std::vector<T> &value){
	  BOOST_FOREACH(const boost::property_tree::ptree::value_type& child,_jsonData.get_child(eleName)){
		value.push_back(child.second.get<T>(""));
	  }
	}

	template<class T>
	void actualRead(const std::string eleName, std::vector<std::vector<T> > &value){
	  BOOST_FOREACH(const boost::property_tree::ptree::value_type& child,_jsonData.get_child(eleName)){
		value.push_back( std::vector<T>() );
		BOOST_FOREACH(const boost::property_tree::ptree::value_type& child_child, child.second){
		  value[value.size()-1].push_back(child_child.second.get<T>(""));
		}
	  }
	}

	template<class T>
	void actualRead(const std::string eleName, std::set<T> &value){
	  BOOST_FOREACH(const boost::property_tree::ptree::value_type& child,_jsonData.get_child(eleName)){
		value.insert(child.second.get<T>(""));
	  }
	}

  private:
	std::ifstream _file;
	std::string _initFilename; 
	std::string _initFileDir;
	bool _printData;
	ptree _jsonData;
  };

}//end of namespace

#endif /* _JSONFILEPASER_H_ */

