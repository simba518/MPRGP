#ifndef _VTKWRITER_H_
#define _VTKWRITER_H_

#include <vector>
#include <string>
#include <fstream>
#include <assertext.h>
#include <Log.h>
using namespace std;

namespace UTILITY{
  
  template<typename T,int N>
  struct VectorXdTmp{
  public:	
	int size()const{
	  return N;
	}
	T &operator [](int i){
	  assert_in(i,0,N-1);
	  return d[i];
	}
	const T &operator[](int i)const{
	  assert_in(i,0,N-1);
	  return d[i];
	}
  private:
    T d[N];
  };

  class endianness{
  public:
	static bool isLittleEndian(){
	  union u { 
		unsigned long l; 
		unsigned char c[sizeof(unsigned long)]; 
	  };
	  u dummy;
	  dummy.l = 1;
	  return dummy.c[0] == 1;
	}
	static void swap2Bytes(unsigned char* &ptr){ 
	  unsigned char tmp;  
	  tmp = ptr[0]; ptr[0] = ptr[1]; ptr[1] = tmp;
	}
	static void swap4Bytes(unsigned char* &ptr){ 
	  unsigned char tmp; 
	  tmp = ptr[0]; ptr[0] = ptr[3]; ptr[3] = tmp;
	  tmp = ptr[1]; ptr[1] = ptr[2]; ptr[2] = tmp; 
	}
	static void swap8Bytes(unsigned char* &ptr){ 
	  unsigned char tmp;
	  tmp = ptr[0]; ptr[0] = ptr[7]; ptr[7] = tmp;
	  tmp = ptr[1]; ptr[1] = ptr[6]; ptr[6] = tmp;
	  tmp = ptr[2]; ptr[2] = ptr[5]; ptr[5] = tmp;
	  tmp = ptr[3]; ptr[3] = ptr[4]; ptr[4] = tmp; 
	}
  };
  template<class T2>
  inline void byteSwap(T2& val){
	size_t n=sizeof(T2);
	unsigned char *p=(unsigned char*)&val;
	switch( n ) 
	  {case 1:return;
	  case 2:endianness::swap2Bytes(p);break;
	  case 4:endianness::swap4Bytes(p);break;
	  case 8:endianness::swap8Bytes(p);break;
	  default:break;}
  }

  template <typename T2>
  inline void vtkWriteBinary(ostream& oss,T2 val){
	if(endianness::isLittleEndian())
	  byteSwap(val);
	oss.write((const char*)&val,sizeof(T2));
  }

  /**
   * @class VTKWriter auxiliary class for helping to write VTK files.
   * http://dunne.uni-hd.de/VisuSimple/documents/vtkfileformat.html
   */
  class VTKWriter{
	
  public:
	VTKWriter(const bool binary=true):
	  _binary(binary),
	  _points(binary ? ios_base::binary : ios_base::out),
	  _cells(binary ? ios_base::binary : ios_base::out),
	  _cellTypes(binary ? ios_base::binary : ios_base::out){
	  
	  clear();
	}

	template<typename VECTOR>
	void addPoints(const VECTOR &points){
	  vector<VectorXdTmp<double,3> > vv;
	  convetVec2VV3d(points,vv);
	  addPoints(vv);
	}

	template<typename VECTOR_3D,typename ALLOCATOR>
	void addPoints(const vector<VECTOR_3D,ALLOCATOR> &points){

	  if (_binary){
		for (int i = 0; i < points.size(); ++i){
		  vtkWriteBinary(_points,points[i][0]);
		  vtkWriteBinary(_points,points[i][1]);
		  vtkWriteBinary(_points,points[i][2]);
		}
	  }else{
		for (int i = 0; i < points.size(); ++i)
		  _points << points[i][0] << " " << points[i][1]<< " " << points[i][2] << endl;
	  }
	  _nPoint += points.size();
	  if (points.size() > 0){
		_dataType = "double";
		if (sizeof(points[0][0])==sizeof(float)){
		  _dataType = "float";
		}else if(sizeof(points[0][0])==sizeof(int)){
		  _dataType = "int";
		}else if(sizeof(points[0][0])!=sizeof(double)){
		  ERROR_LOG("data type for the point is not support!");
		}
	  }
	}

	template<typename VECTOR>
	void addLines(const VECTOR &lines){
	  vector<VectorXdTmp<int,2> > vv;
	  convetVec2VV2d(lines,vv);
	  addLines(vv);
	}

	template<typename VECTOR_2I,typename ALLOCATOR>
	void addLines(const vector<VECTOR_2I,ALLOCATOR> &lines){

	  if (lines.size() > 0)
		assert_eq(lines[0].size(),2);
	  if (_binary){
		for (int i = 0; i < lines.size(); ++i){
		  vtkWriteBinary(_cells,2);
		  vtkWriteBinary(_cells,lines[i][0]);
		  vtkWriteBinary(_cells,lines[i][1]);
		  vtkWriteBinary(_cellTypes,3);
		}
	  }else{
		for (int i = 0; i < lines.size(); ++i){
		  _cells<< 2 <<" "<< lines[i][0] << " " << lines[i][1]<< endl;
		  _cellTypes << 3 << " ";
		}
	  }
	  _nCells += lines.size();
	  _nCellIndexes += lines.size()*(2+1);
	}

	template<typename VECTOR>
	void addTriangles(const VECTOR &faces){
	  vector<VectorXdTmp<int,3> > vv;
	  convetVec2VV3d(faces,vv);
	  addTriangles(vv);
	}

	template<typename VECTOR_3I,typename ALLOCATOR>
	void addTriangles(const vector<VECTOR_3I,ALLOCATOR> &faces){

	  if (faces.size() > 0)
		assert_eq(faces[0].size(),3);
	  if (_binary){
		for (int i = 0; i < faces.size(); ++i){
		  vtkWriteBinary(_cells,3);
		  vtkWriteBinary(_cells,faces[i][0]);
		  vtkWriteBinary(_cells,faces[i][1]);
		  vtkWriteBinary(_cells,faces[i][2]);
		  vtkWriteBinary(_cellTypes,5);
		}
	  }else{
		for (int i = 0; i < faces.size(); ++i){
		  _cells<< 3 <<" "<< faces[i][0] << " " << faces[i][1]<< " " << faces[i][2] << endl;
		  _cellTypes << 5 << " ";
		}
	  }
	  _nCells += faces.size();
	  _nCellIndexes += faces.size()*(3+1);
	}

	template<typename VECTOR>
	void addTets(const VECTOR &tets){
	  vector<VectorXdTmp<int,4> > vv;
	  convetVec2VV3d(tets,vv);
	  addTets(vv);
	}

	template<typename VECTOR_4I,typename ALLOCATOR>
	void addTets(const vector<VECTOR_4I,ALLOCATOR> &tets){

	  if (tets.size() > 0)
		assert_eq(tets[0].size(),4);

	  if (_binary){
		for (int i = 0; i < tets.size(); ++i){
		  vtkWriteBinary(_cells,4);
		  vtkWriteBinary(_cells,tets[i][0]);
		  vtkWriteBinary(_cells,tets[i][1]);
		  vtkWriteBinary(_cells,tets[i][2]);
		  vtkWriteBinary(_cells,tets[i][3]);
		  vtkWriteBinary(_cellTypes,10);
		}
	  }else{
		for (int i = 0; i < tets.size(); ++i){
		  _cells<<4<<" "<<tets[i][0]<<" "<<tets[i][1]<<" "<<tets[i][2]<<" "<<tets[i][3]<<endl;
		  _cellTypes << 10 << " ";
		}
	  }

	  _nCells += tets.size();
	  _nCellIndexes += tets.size()*(4+1);
	}

	bool write(const string filename){

	  ofstream outf(filename.c_str(),_binary?ios_base::binary:ios_base::out);
	  bool succ = outf.is_open();
	  if (succ){
		if (_nCells <= 0)
		  createVertexCells();
		outf << head() << endl;
		outf << "\nPOINTS " << _nPoint<<" "<< _dataType << endl;
		outf << _points.str() << endl;
		outf << "\nCELLS " << _nCells<<" " << _nCellIndexes << endl;
		outf << _cells.str() << endl;
		outf << "CELL_TYPES " << _nCells << endl;
		outf << _cellTypes.str() << endl;
	  }
	  outf.close();
	  return succ;
	}

	void clear(){
	  _nPoint = 0;
	  _nCells = 0;
	  _nCellIndexes = 0;
	  _points.clear();
	  _cellTypes.clear();
	  _cells.clear();
	}
	
  protected:
	string head()const{
	  const string s1="# vtk DataFile Version 3.1\n";
	  const string s2="generated by UTILITY::VTKWriter in Utility library,author: simba\n";
	  const string s3=(_binary ? "BINARY\n":"ASCII\n");
	  const string s4="DATASET UNSTRUCTURED_GRID\n";
	  return s1+s2+s3+s4;
	}

	template<typename VECTOR,typename VECTOR3D>
	void convetVec2VV3d(const VECTOR &vec,VECTOR3D &vvec)const{

	  assert_eq(vec.size()%3,0);
	  vvec.resize(vec.size()/3);
	  for (size_t i = 0; i < vvec.size(); ++i){
		vvec[i][0] = vec[i*3+0];
		vvec[i][1] = vec[i*3+1];
		vvec[i][2] = vec[i*3+2];
	  }
	}

	template<typename VECTOR,typename VECTOR2D>
	void convetVec2VV2d(const VECTOR &vec,VECTOR2D &vvec)const{

	  assert_eq(vec.size()%2,0);
	  vvec.resize(vec.size()/2);
	  for (size_t i = 0; i < vvec.size(); ++i){
		vvec[i][0] = vec[i*2+0];
		vvec[i][1] = vec[i*2+1];
	  }
	}

	void createVertexCells(){
	  if (_binary){
		for (int i = 0; i < _nPoint; ++i){
		  vtkWriteBinary(_cells,1);
		  vtkWriteBinary(_cells,i);
		  vtkWriteBinary(_cellTypes,1);
		}
	  }else{
		for (int i = 0; i < _nPoint; ++i){
		  _cells<<1<< " " << i << endl;
		  _cellTypes << 1 << " ";
		}
	  }

	  _nCells += _nPoint;
	  _nCellIndexes += _nPoint*2;
	}
	
  private:
	const bool _binary;
	size_t _nPoint,_nCells,_nCellIndexes;
	ostringstream _points,_cells,_cellTypes;
	string _dataType;
  };

  // other methods
  template <typename OS, typename INT>
  void tet2vtk(OS &os,const double *node, size_t node_num,const INT *tet, size_t tet_num){

	os << "# vtk DataFile Version 2.0\nTET\nASCII\n\nDATASET UNSTRUCTURED_GRID\n";
	os << "POINTS " << node_num << " double\n";
	for(size_t i = 0; i < node_num; ++i)
	  os << node[i*3+0] << " " << node[i*3+1] << " " << node[i*3+2] << "\n";

	os << "CELLS " << tet_num << " " << tet_num*5 << "\n";
	for(size_t i = 0; i < tet_num; ++i)
	  os << 4 << "  "
		 << tet[i*4+0] << " " << tet[i*4+1] << " "
		 << tet[i*4+2] << " " << tet[i*4+3] << "\n";
	os << "CELL_TYPES " << tet_num << "\n";
	for(size_t i = 0; i < tet_num; ++i)
	  os << 10 << "\n";
  }

  template <typename OS, typename Iterator, typename INT>
  void vtk_data(OS &os, Iterator first, INT size, const char *value_name, const char *table_name = "my_table"){
	os << "SCALARS " << value_name << " double\nLOOKUP_TABLE " << table_name << "\n";
	for(size_t i = 0; i < size; ++i, ++first)
	  os << *first << "\n";
  }

  template <typename OS, typename Iterator, typename INT>
  void cell_data(OS &os, Iterator first, INT size, const char *value_name, const char *table_name = "my_table"){
	os << "CELL_DATA " << size << "\n";
	vtk_data(os, first, size, value_name, table_name);
  }
  
}//end of namespace

#endif /*_VTKWRITER_H_*/
