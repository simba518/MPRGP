#ifndef _ASSERTMACROSEXT_H_
#define _ASSERTMACROSEXT_H_

/**
 * This is some extentions for the assertion macros of the CPPUNIT_TEST
 * framework, which is more powerfull and more easy to use.
 */

#include <sstream>
#include <iomanip>
#include <cppunit/TestAssert.h>

// assert "value_a == value_b", if not, print the values and abort.
#define ASSERT_EQ(value_a, value_b)										\
  if (!(value_a == value_b))											\
	{																	\
	  std::ostringstream stringStream;									\
	  stringStream << __STRING(value_a) <<" = " << value_a<< std::endl; \
	  stringStream << "- ";												\
	  stringStream <<__STRING(value_b) <<" = " << value_b<< std::endl;	\
	  CPPUNIT_ASSERT_MESSAGE(stringStream.str(),value_a == value_b);	\
	}														

// assert "||value_a - value_b||<=tol", if not, print the values and abort.
#define ASSERT_EQ_TOL(value_a, value_b,tol)								\
  if (std::abs(value_a - value_b)>tol)										\
	{																	\
	  std::ostringstream stringStream;									\
	  stringStream << __STRING(value_a) <<" = " << value_a<< std::endl; \
	  stringStream << "- ";												\
	  stringStream <<__STRING(value_b) <<" = " << value_b<< std::endl;	\
	  CPPUNIT_ASSERT_MESSAGE(stringStream.str(),std::abs(value_a-value_b)<=tol); \
	}														

// assert "value_a != value_b", if not, print the values and abort.
#define ASSERT_NE(value_a, value_b)										\
  if (value_a == value_b)												\
	{																	\
	  std::ostringstream stringStream;									\
	  stringStream << __STRING(value_a) <<" = " << value_a<< std::endl; \
	  stringStream << "- ";												\
	  stringStream <<__STRING(value_b) <<" = " << value_b<< std::endl;	\
	  CPPUNIT_ASSERT_MESSAGE(stringStream.str(),value_a != value_b);	\
	}														

// assert "value_a >= value_b", if not, print the values and abort.
#define ASSERT_GE(value_a, value_b)										\
  if (value_a < value_b)												\
	{																	\
	  std::ostringstream stringStream;									\
	  stringStream << __STRING(value_a) <<" = " << value_a<< std::endl; \
	  stringStream << "- ";												\
	  stringStream <<__STRING(value_b) <<" = " << value_b<< std::endl;	\
	  CPPUNIT_ASSERT_MESSAGE(stringStream.str(),value_a >= value_b);	\
	}														

// assert "value_a > value_b", if not, print the values and abort.
#define ASSERT_GT(value_a, value_b)										\
  if (value_a <= value_b)												\
	{																	\
	  std::ostringstream stringStream;									\
	  stringStream << __STRING(value_a) <<" = " << value_a<< std::endl; \
	  stringStream << "- ";												\
	  stringStream <<__STRING(value_b) <<" = " << value_b<< std::endl;	\
	  CPPUNIT_ASSERT_MESSAGE(stringStream.str(),value_a > value_b);		\
	}														

// assert "value_a <= value_b", if not, print the values and abort.
#define ASSERT_LE(value_a, value_b)										\
  if (value_a > value_b)												\
	{																	\
	  std::ostringstream stringStream;									\
	  stringStream << __STRING(value_a) <<" = " << value_a<< std::endl; \
	  stringStream << "- ";												\
	  stringStream <<__STRING(value_b) <<" = " << value_b<< std::endl;	\
	  CPPUNIT_ASSERT_MESSAGE(stringStream.str(),value_a <= value_b);	\
	}														

// assert "value_a < value_b", if not, print the values and abort.
#define ASSERT_LT(value_a, value_b)										\
  if ( !(value_a < value_b) )												\
	{																	\
	  std::ostringstream stringStream;									\
	  stringStream << __STRING(value_a) <<" = " << value_a<< std::endl; \
	  stringStream << "- ";												\
	  stringStream <<__STRING(value_b) <<" = " << value_b<< std::endl;	\
	  CPPUNIT_ASSERT_MESSAGE(stringStream.str(),value_a < value_b);		\
	}														

// @brief assert "vector_a == vector_b", if not, abort. 
// @note the type of vector_a and vector_a could be int*, float*, double* or
// std::vector. 
// @note len is the length of the vectors. and both vector_a and vector_b should
// have the same length, which won't check here.
// @note the vector here should not be too large, as we want to print them when
// the assertion failed.
#define ASSERT_EQ_SMALL_VEC(vector_a,vector_b,len)						\
  for (int i = 0; i < len; ++i)											\
	{																	\
	  if (vector_a[i] != vector_b[i])									\
		{																\
		  std::ostringstream v_a;										\
		  std::ostringstream v_b;										\
		  std::ostringstream v_a_b;										\
		  v_a_b << "- |" << __STRING(vector_a) <<"["<<i<<"] - " << __STRING(vector_b) <<"["<<i<<"]|=" << std::abs(vector_a[i]-vector_b[i]); \
		  v_a << "- " << __STRING(vector_a) << " = ( ";					\
		  v_b << "- " << __STRING(vector_b) << " = ( ";					\
		  for (int j = 0; j < len; ++j)									\
			{															\
			  v_a << std::setprecision(5) << vector_a[j]<<" ";			\
			  v_b << std::setprecision(5) << vector_b[j]<<" ";			\
			}															\
		  v_a << " ) " << endl;											\
		  v_b << " ) " << endl;											\
		  v_a << v_b.str() << endl << v_a_b.str() << endl;				\
		  CPPUNIT_ASSERT_MESSAGE(v_a.str(),vector_a == vector_b);		\
		  break;														\
		}																\
	}																	\


#define ASSERT_EQ_SMALL_VEC_TOL(vector_a,vector_b,len, tol)				\
  for (int i = 0; i < len; ++i)											\
	{																	\
	  if (std::abs(vector_a[i]-vector_b[i]) > tol)							\
		{																\
		  std::ostringstream v_a;										\
		  std::ostringstream v_b;										\
		  std::ostringstream v_a_b;										\
		  v_a_b << "- |" << __STRING(vector_a) <<"["<<i<<"] - " << __STRING(vector_b) <<"["<<i<<"]|=" << std::abs(vector_a[i]-vector_b[i]); \
		  v_a << "- " << __STRING(vector_a) << " = ( ";					\
		  v_b << "- " << __STRING(vector_b) << " = ( ";					\
		  for (int j = 0; j < len; ++j)									\
			{															\
			  v_a << std::setprecision(5) << vector_a[j]<<" ";			\
			  v_b << std::setprecision(5) << vector_b[j]<<" ";			\
			}															\
		  v_a << " ) " << endl;											\
		  v_b << " ) " << endl;											\
		  v_a << v_b.str() << endl;										\
		  v_a << v_a_b.str() << endl;									\
		  CPPUNIT_ASSERT_MESSAGE(v_a.str(),vector_a == vector_b);		\
		  break;														\
		}																\
	}																	\

// @brief assert "mat_a == mat_b", if not, abort. 
// @note mat_a and mat_b usually should be MatrixXd.
// @note the mat here should not be too large, as we want to print them when
// the assertion failed.
#define ASSERT_EQ_SMALL_MAT(mat_a,mat_b)								\
  if (mat_a.rows()!=mat_b.rows() || mat_a.cols()!=mat_b.cols()){		\
																		\
	std::ostringstream msg;												\
	msg<<__STRING(mat_a) << ".dimension() != ";							\
	msg<<__STRING(mat_b) << ".dimension()\n-";							\
	msg<<__STRING(mat_a)<<":("<<mat_a.rows()<<","<<mat_a.cols()<<")\n-"; \
	msg<<__STRING(mat_b)<<":("<<mat_b.rows()<<","<<mat_b.cols()<<")\n";	\
	CPPUNIT_ASSERT_MESSAGE(msg.str(),false);							\
																		\
  }else{																\
  																		\
	for (int i = 0; i < mat_a.rows(); ++i)								\
	  for (int j = 0; j < mat_a.cols(); ++j)							\
		{																\
		  if (mat_a(i,j) != mat_b(i,j))									\
			{															\
			  std::ostringstream msg;									\
		  std::ostringstream v_a_b;										\
			  msg << __STRING(mat_b)<<":\n" << mat_b << endl << endl;	\
			  msg << __STRING(mat_a)<<":\n" << mat_a << endl << endl;	\
			  msg << "- |" << __STRING(mat_a) <<"("<<i<<","<<j<<") - " << __STRING(mat_b) <<"("<<i<<","<<j<<")|=" << std::abs(mat_a(i,j)-mat_b(i,j)) << endl; \
			  CPPUNIT_ASSERT_MESSAGE(msg.str(),mat_a == mat_b);			\
			  break;													\
			}															\
		}																\
   }

#define ASSERT_EQ_SMALL_MAT_TOL(mat_a,mat_b,tol)						\
  if (mat_a.rows()!=mat_b.rows() || mat_a.cols()!=mat_b.cols()){		\
																		\
	std::ostringstream msg;												\
	msg<<__STRING(mat_a) << ".dimension() != ";							\
	msg<<__STRING(mat_b) << ".dimension()\n-";							\
	msg<<__STRING(mat_a)<<":("<<mat_a.rows()<<","<<mat_a.cols()<<")\n-"; \
	msg<<__STRING(mat_b)<<":("<<mat_b.rows()<<","<<mat_b.cols()<<")\n";	\
	CPPUNIT_ASSERT_MESSAGE(msg.str(),false);							\
																		\
  }else{																\
  																		\
	for (int i = 0; i < mat_a.rows(); ++i)								\
	  for (int j = 0; j < mat_a.cols(); ++j)							\
		{																\
		  if (std::abs(mat_a(i,j)-mat_b(i,j)) > tol)						\
			{															\
			  std::ostringstream msg;									\
			  msg << __STRING(mat_a)<<":\n" << mat_a << endl << endl;	\
			  msg << __STRING(mat_b)<<":\n" << mat_b << endl << endl;	\
			  msg << "- |" << __STRING(mat_a) <<"("<<i<<","<<j<<") - " << __STRING(mat_b) <<"("<<i<<","<<j<<")|=" << std::abs(mat_a(i,j)-mat_b(i,j)) << endl; \
			  CPPUNIT_ASSERT_MESSAGE(msg.str(),mat_a == mat_b);			\
			  break;													\
			}															\
		}																\
   }

#endif /* _ASSERTMACROSEXT_H_ */
