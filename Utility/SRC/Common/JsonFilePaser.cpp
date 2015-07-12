#include "JsonFilePaser.h"

#ifdef WIN32
#define _interlockedbittestandreset _interlockedbittestandreset_NAME_CHANGED_TO_AVOID_MSVS2008_ERROR
#define _interlockedbittestandset _interlockedbittestandset_NAME_CHANGED_TO_AVOID_MSVS2008_ERROR
#endif

#include <boost/property_tree/json_parser.hpp>
using namespace UTILITY;

JsonFilePaser::JsonFilePaser(bool print){
  _printData = print;
}

JsonFilePaser::~JsonFilePaser(){
  close();
}

bool JsonFilePaser::open(const std::string filename){

  bool succ = false;
  close();
  _initFilename = filename;
  _file.open(filename.c_str());
  ERROR_LOG_COND("failed to open init file: "<<filename,_file.is_open());
  if(_file.is_open()){
	try{
	  json_parser::read_json(_file,_jsonData);
	  succ = read("init_file_dir",_initFileDir);
	  _initFileDir = replaceHomeDir(_initFileDir);
	}catch(std::exception& ex){
	  ERROR_LOG("file: "<< filename << ex.what());
	}
  }
  return succ;
}
void JsonFilePaser::close(){ if(_file.is_open()) _file.close();}

std::string JsonFilePaser::getFileDir()const{
  return _initFileDir;
}
std::string JsonFilePaser::getFileName()const{
  return _initFilename;
}

bool JsonFilePaser::readFilePath(const std::string eleName, std::string &filePath,const bool checkFileExist){
  bool succ = false;
  if( read(eleName,filePath) ){
	filePath = getFileDir()+filePath;
	succ = fileIsExist(eleName,filePath,checkFileExist);
  }
  return succ;
}
bool JsonFilePaser::readFilePath(const std::string eleName, std::vector<std::string> &filePathes,const bool checkFileExist){

  bool succ = false;
  if(read<std::vector<std::string> >(eleName,filePathes)){
	succ = true;
	for (size_t i = 0; i < filePathes.size(); ++i){
	  filePathes[i] = getFileDir()+filePathes[i];
	  succ &= fileIsExist(eleName,filePathes[i],checkFileExist);
	}
  }
  return succ;
}
std::string JsonFilePaser::replaceHomeDir(const std::string &path,const bool checkPathExist){
  /// replace the first "~" as the home directory
  std::string newPath = path;
  if(path.size() > 0 && path[0] == '~')
	newPath = getenv("HOME") + newPath.erase(0,1);
  dirIsExist("",newPath,checkPathExist);
  return newPath;
}
	
bool JsonFilePaser::fileIsExist(const std::string eleName,const std::string filePath,bool check)const{
  bool exist = true;
  if(check){
	exist = (boost::filesystem::exists(filePath) && 
			 !(boost::filesystem::is_directory(filePath)));
	WARN_LOG_COND("file '"<< filePath <<"' is not existed! (in node '" << eleName <<"' )" ,exist);
  }
  return exist;
}
bool JsonFilePaser::dirIsExist(const std::string eleName,const std::string dirPath,bool check)const{

  bool exist = true;
  if(check){
		boost::filesystem::path mypath(dirPath.c_str());
		exist = boost::filesystem::is_directory(mypath);
		WARN_LOG_COND("dir '"<< dirPath <<"' is not existed! (in node '" << eleName <<"' )" ,exist);
	}
  return exist;
}
