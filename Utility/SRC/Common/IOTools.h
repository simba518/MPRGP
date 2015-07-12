#ifndef _IOTOOLS_H_
#define _IOTOOLS_H_

#include <fstream>
#include <vector>
#include <set>
#include <assertext.h>
#include <Log.h>

namespace UTILITY{

  enum IO_TYPE{BINARY,TEXT};

  template<class T>
  inline bool loadText(std::ifstream &in,T *data,const int len,const std::string fname=""){

	assert(data);
	assert_ge(len,0);
	bool succ = false;
	if (in.is_open()){
	  succ = true;
	  for (int i = 0; i < len; ++i){
		if ( in.eof() || !(in >> (data[i])) ){
		  ERROR_LOG("failed to read data,"<<fname<<",i="<<i);
		  succ = false;
		  break;
		}
	  }
	}
	return succ;
  }

  template<class T>
  inline bool writeText(std::ofstream &out,const T *data,const int len,const std::string space="\t",const std::string fname=""){

	bool succ = false;
	if (out.is_open() && data != NULL){
	  succ = true;
	  for ( int i = 0; i < len ; i++){
		if ( !(out <<data[i]<<space) ) {
		  ERROR_LOG("failed to save data,"<<fname<<",i="<<i);
		  succ = false;
		  break;
		}
	  }
	}
	return succ;
  }

  template<class T>
  inline bool loadBinary(std::ifstream &in,T *data,const int len,const std::string fname=""){
		
	bool succ = true;
	const int nRead = in.read((char*)(data), sizeof(T)*len).gcount();
	if (nRead != len*(int)sizeof(T)){
	  ERROR_LOG("failed to read data,"<<fname);
	  succ = false;
	}
	return succ;
  }

  template<class T>
  inline bool writeBinary(std::ofstream &out,const T *data,const int len,const std::string fname=""){

	bool succ = true;
	out.write((char*)(data),sizeof(T)*len);
	if (out.fail()){
	  ERROR_LOG("failed to save data,"<<fname);
	  succ = false;
	}
	return succ;
  }

  template<class T>
  inline bool load(std::ifstream &in,T *data,const int len,const IO_TYPE io_type,const std::string fname=""){
	if(BINARY == io_type){
	  return loadBinary(in,data,len,fname);
	}else{
	  return loadText(in,data,len,fname);
	}
  }

  template<class T>
  inline bool write(std::ofstream &out,const T *data,const int len,const IO_TYPE io_type,const std::string space="\t",const std::string fname=""){
	if(BINARY == io_type){
	  return writeBinary(out,data,len,fname);
	}else{
	  return writeText(out,data,len,space,fname);
	}
  }

  inline bool openInputFile(std::ifstream &in,const std::string fname,IO_TYPE io_type){

	if(BINARY == io_type){
	  in.open(fname.c_str(),std::ios_base::in|std::ios_base::binary);
	}else{
	  in.open(fname.c_str());
	}
	const bool succ = in.is_open();
	ERROR_LOG_COND("failed to open file: "<<fname,succ);
	return succ;
  }

  inline bool openOutputFile(std::ofstream &out,const std::string fname,IO_TYPE io_type){

	if(BINARY == io_type){
	  out.open(fname.c_str(),std::ios_base::out|std::ios_base::binary);
	}else{
	  out.open(fname.c_str());
	}
	const bool succ = out.is_open();
	ERROR_LOG_COND("failed to open file: "<<fname,succ);
	return succ;
  }

  template<class VECTOR>
  inline bool loadVec(const std::string fname,VECTOR &data,IO_TYPE io_type=BINARY){

	bool succ = false;
	std::ifstream in;	
	if( openInputFile(in,fname,io_type) ){
	  int len = 0;
	  if(load(in,&len,1,io_type,fname)){
		if (len <= 0){
		  ERROR_LOG("data length is zero in file: "<<fname);
		  return false;
		}
		data.resize(len);
		succ = load(in,&data[0],len,io_type,fname);
	  }
	}
	in.close();
	return succ;
  }

  template<class T>
  inline bool load(const std::string fname,std::set<T> &data,IO_TYPE io_type=BINARY){
	std::vector<T> dv;
	const bool succ = loadVec(fname,dv,io_type);
	data.clear();
	if (succ)
	  for (size_t i = 0; i < dv.size(); ++i)
		data.insert(dv[i]);
	return succ;
  }

  template<class VECTOR>
  inline bool writeVec(const std::string fname,const VECTOR &data,IO_TYPE io_type=BINARY,const std::string space="\t"){
	
	bool succ = false;
	std::ofstream out;
	if( openOutputFile(out,fname,io_type) ){
	  const int len = (int)data.size();
	  if(write(out,&len,1,io_type,space,fname)){
		if (len <= 0){
		  ERROR_LOG("data length is zero for writing file: "<<fname);
		  return false;
		}
		succ = write(out,&data[0],len,io_type,space,fname);
	  }
	}
	return succ;
  }

  template<class MATRIX>
  inline bool loadMat(const std::string fname,MATRIX &data,IO_TYPE io_type=BINARY){

	bool succ = false;
	std::ifstream in;	
	if( openInputFile(in,fname,io_type) ){
	  int rows,cols;
	  if(load(in,&rows,1,io_type,fname) && load(in,&cols,1,io_type,fname) ){
		if (rows <=0 || cols <= 0){
		  ERROR_LOG("data dimension in file: "<<fname << " is less than zero "<<"rows="<<rows<<",cols="<<cols);
		  return false;
		}
		data.resize(rows,cols);
		succ = load(in,&data(0,0),rows*cols,io_type,fname);
	  }
	}
	in.close();
	return succ;
  }

  template<class MATRIX>
  inline bool writeMat(const std::string fname,const MATRIX &data,IO_TYPE io_type=BINARY,const std::string space="\t"){
	
	bool succ = false;
	std::ofstream out;
	if( openOutputFile(out,fname,io_type) ){
	  const int rows = (int)data.rows();
	  const int cols = (int)data.cols();
	  if(write(out,&rows,1,io_type,fname,space) && write(out,&cols,1,io_type,fname,space)){
		if (rows <=0 || cols <= 0){
		  ERROR_LOG("data dimension is less than zero "<<"rows="<<rows<<",cols="<<cols);
		  return false;
		}
		succ = write(out,&data(0,0),rows*cols,io_type,space,fname);
	  }
	}
	return succ;
  }
}

#endif /* _IOTOOLS_H_ */
