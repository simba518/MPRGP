#ifndef _MATRIXIO_H_
#define _MATRIXIO_H_

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <IOTools.h>

namespace EIGEN3EXT{

  template<class T>
  inline bool load(const std::string fname,Eigen::Matrix<T,1,-1> &data,UTILITY::IO_TYPE io_type=UTILITY::BINARY){
	return UTILITY::loadVec(fname,data,io_type);
  }

  template<class T>
  inline bool load(const std::string fname,Eigen::Matrix<T,-1,1> &data,UTILITY::IO_TYPE io_type=UTILITY::BINARY){
	return UTILITY::loadVec(fname,data,io_type);
  }

  template<class T>
  inline bool write(const std::string fname,const Eigen::Matrix<T,1,-1> &data,UTILITY::IO_TYPE io_type=UTILITY::BINARY,const std::string space="\t"){
	return UTILITY::writeVec(fname,data,io_type);
  }

  template<class T>
  inline bool write(const std::string fname,const Eigen::Matrix<T,-1,1> &data,UTILITY::IO_TYPE io_type=UTILITY::BINARY,const std::string space="\t"){
	return UTILITY::writeVec(fname,data,io_type);
  }

  template<class T>
  inline bool load(const std::string fname,Eigen::Matrix<T,-1,-1> &data,UTILITY::IO_TYPE io_type=UTILITY::BINARY){
    return UTILITY::loadMat(fname,data,io_type);
  }

  template<class T>
  inline bool write(const std::string fname,const Eigen::Matrix<T,-1,-1> &data,UTILITY::IO_TYPE io_type=UTILITY::BINARY,const std::string space="\t"){
	return UTILITY::writeMat(fname,data,io_type,space);
  }

  template<class T>
  inline bool load(const std::string fname,std::vector<Eigen::Matrix<T,-1,1> > &data,UTILITY::IO_TYPE io_type=UTILITY::BINARY){

	bool succ = false;
	std::ifstream in;	
	if( openInputFile(in,fname,io_type) ){
	  int Total,subLen;
	  if(load(in,&subLen,1,io_type,fname) && load(in,&Total,1,io_type,fname) ){
		assert_gt(Total,0);
		assert_gt(subLen,0);
		data.resize(Total);
		for (size_t i = 0; i < data.size(); ++i){
		  data[i].resize(subLen);
		  succ = load(in,&(data[i][0]),subLen,io_type,fname);
		  if(!succ){
			break;
		  }
		}
	  }
	}
	in.close();
	return succ;
  }

  template<class T>
  inline bool write(const std::string fname,const std::vector<Eigen::Matrix<T,-1,1> > &data,UTILITY::IO_TYPE io_type=UTILITY::BINARY,const std::string space="\t"){
	
	bool succ = false;
	std::ofstream out;
	if( openOutputFile(out,fname,io_type) ){
	  const int Total = (int)data.size();
	  assert_gt(Total,0);
	  const int subLen = (int)data[0].size();
	  assert_gt(subLen,0);
	  if(write(out,&subLen,1,io_type,fname,space) && write(out,&Total,1,io_type,fname,space)){
		for (int i = 0; i < Total; ++i){
		  succ = write(out,&data[i][0],subLen,io_type,space,fname);
		  if(!succ)
			break;
		}
	  }
	}
	return succ;
  }

  template <typename T> 
  bool loadBinary(int &rows,int &cols,
				  std::vector<int> &row_ids,std::vector<int> &col_ids,
				  std::vector<T> &values,const std::string &f_name){

	
	// open file
	std::ifstream in_f;
	in_f.open(f_name.c_str(),std::ios_base::in|std::ios_base::binary);
	bool succ = in_f.is_open();

	// read data
	if (succ){
	  
	  int nz = 0;
	  in_f.read((char*)(&nz),sizeof(int));    // read number of nonzeros
	  in_f.read((char*)(&rows),sizeof(int));  // read rows
	  in_f.read((char*)(&cols),sizeof(int));  // read cols
	  succ = !(in_f.fail());
	  
	  assert_ge(rows,0);
	  assert_ge(cols,0);
	  assert_ge(rows*cols,nz);
	  assert_ge(nz,0);
	  
	  row_ids.resize(nz);
	  col_ids.resize(nz);
	  values.resize(nz);
	  
	  if (succ && nz > 0){
		in_f.read((char*)(&row_ids[0]), sizeof(int)*(row_ids.size()));
		in_f.read((char*)(&col_ids[0]), sizeof(int)*(col_ids.size()));
		in_f.read((char*)(&values[0]), sizeof(T)*(values.size()));
		succ = !(in_f.fail());
	  }
	  in_f.close(); // close file
	  ERROR_LOG_COND("failed to load the sparse matrix from: "<<f_name,succ);
	}else{
	  ERROR_LOG("failed to open file: "<<f_name);
	}
	return succ; 
  }

  template <typename T> 
  bool loadText(int &rows,int &cols,
				  std::vector<int> &row_ids,std::vector<int> &col_ids,
				  std::vector<T> &values,const std::string &f_name){
	/// @todo
	ERROR_LOG("undefined function is called.");
	return false;
  }

  template <typename T> 
  bool writeBinary(const int rows,const int cols,
				  const std::vector<int> &row_ids,const std::vector<int> &col_ids,
				  const std::vector<T> &values,const std::string &f_name){

	// open file
	std::ofstream out_f;
	out_f.open(f_name.c_str(),std::ios_base::out|std::ios_base::binary);
	bool succ = out_f.is_open();

	// write data
	if (succ){
	  
	  const int nz = values.size();
	  out_f.write((char*)(&nz),sizeof(int));    // write nz
	  out_f.write((char*)(&rows),sizeof(int));  // write rows
	  out_f.write((char*)(&cols),sizeof(int));  // write cols
	  succ = !(out_f.fail());
	  if (succ && nz > 0){

		assert_eq ((int)row_ids.size(), nz);
		assert_eq ((int)col_ids.size(), nz);
		out_f.write((char*)(&row_ids[0]), sizeof(int)*(row_ids.size()));
		out_f.write((char*)(&col_ids[0]), sizeof(int)*(col_ids.size()));
		out_f.write((char*)(&values[0]), sizeof(T)*(values.size()));
		succ = !(out_f.fail());
	  }
	  out_f.close(); // close file
	  ERROR_LOG_COND("failed to save the sparse matrix to: "<<f_name,succ);
	}else{
	  ERROR_LOG("failed to open file: "<<f_name);
	}
	return succ; 
  }

  template <typename T> 
  bool writeText(const int rows,const int cols,
				  const std::vector<int> &row_ids,const std::vector<int> &col_ids,
				  const std::vector<T> &values,const std::string &f_name){
	/// @todo
	ERROR_LOG("undefined function is called.");
	return false;
  }

  template <typename T> 
  bool load(Eigen::SparseMatrix<T>&S,const std::string f_name,UTILITY::IO_TYPE io_type=UTILITY::BINARY){

	// read from file
	bool succ = false;
	std::vector<int> row_ids, col_ids;
	std::vector<T> values;
	int rows, cols;
	if (io_type == UTILITY::TEXT){
	  succ = loadText(rows,cols,row_ids,col_ids,values,f_name);
	}else{
	  succ = loadBinary(rows,cols,row_ids,col_ids,values,f_name);
	}

	// initialize the sparse matrix
	if (succ){

	  assert_ge(rows,0);
	  assert_ge(cols,0);
	  assert_eq(row_ids.size(),values.size());
	  assert_eq(col_ids.size(),values.size());

	  std::vector<Eigen::Triplet<T> > triplets;
	  triplets.reserve(values.size());
	  for (int i = 0; i < (int)values.size(); ++i){
		assert_in(row_ids[i],0,rows-1);
		assert_in(col_ids[i],0,cols-1);
		triplets.push_back( Eigen::Triplet<T>(row_ids[i],col_ids[i],values[i]) );
	  }
	  S.resize(rows,cols);
	  S.reserve(values.size());
	  S.setFromTriplets(triplets.begin(), triplets.end());
	}
	return succ;
  }

  template <typename T> 
  bool write(const Eigen::SparseMatrix<T>&S,const std::string f_name,UTILITY::IO_TYPE io_type=UTILITY::BINARY){

	bool succ = false;
	// the rows, cols and values
	const int nz = S.nonZeros();
	std::vector<int> row_ids(nz), col_ids(nz);
	std::vector<T> values(nz);
	int i = 0;
	for (int k=0; k<S.outerSize(); ++k){ 
	  for (typename Eigen::SparseMatrix<T>::InnerIterator it(S,k); it; ++it){
		row_ids[i] = it.row();
		col_ids[i] = it.col();
		values[i] = it.value();
		i ++;
	  }
	}

	// save to file
	const int rows = S.rows();
	const int cols = S.cols();
	if (io_type == UTILITY::TEXT){
	  succ = writeText(rows,cols,row_ids,col_ids,values,f_name);
	}else{
	  succ = writeBinary(rows,cols,row_ids,col_ids,values,f_name);
	}
	return succ;
  }

}

#endif /* _MATRIXIO_H_ */
