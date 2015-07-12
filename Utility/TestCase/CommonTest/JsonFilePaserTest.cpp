#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <JsonFilePaser.h>
#include <iostream>

struct JsonFilePaserTestInit{
  JsonFilePaserTestInit(){
	jsonFilePath = std::string(TEST_DATA_DIR)+"json.in";
	TEST_ASSERT_MSG("failed to open file:"<<jsonFilePath,jsonFile.open(jsonFilePath));
  }
  std::string jsonFilePath;
  UTILITY::JsonFilePaser jsonFile;
};

BOOST_AUTO_TEST_SUITE(JsonFilePaser)

BOOST_FIXTURE_TEST_CASE(readFloat, JsonFilePaserTestInit){
  float f;
  TEST_ASSERT(jsonFile.read("floatnumber",f));
  ASSERT_EQ(f,100.0);
}

BOOST_FIXTURE_TEST_CASE(readVector, JsonFilePaserTestInit){

  std::vector<int> v, cv;
  cv.push_back(1);
  cv.push_back(2);
  cv.push_back(3);
  TEST_ASSERT(jsonFile.read("intArray",v));
  ASSERT_EQ(v.size(),cv.size());
  ASSERT_EQ_SMALL_VEC(v,cv,v.size());
}

BOOST_FIXTURE_TEST_CASE(readStrVector, JsonFilePaserTestInit){

  std::vector<std::string> v;
  TEST_ASSERT(jsonFile.read("filenames",v));
  ASSERT_EQ(v.size(),3);
  ASSERT_EQ(v[0], "f0.txt");
  ASSERT_EQ(v[1], "f1.txt");
  ASSERT_EQ(v[2], "f2.txt");
}

BOOST_FIXTURE_TEST_CASE(readFilePath, JsonFilePaserTestInit){

  std::string filename;
  TEST_ASSERT(jsonFile.readFilePath("filename",filename,false));
  ASSERT_NE(filename,std::string("f.txt"));
  ASSERT_EQ(filename,jsonFile.getFileDir()+std::string("f.txt"));

  std::vector<std::string> filenames;
  TEST_ASSERT(jsonFile.readFilePath("filenames",filenames,false));
  ASSERT_EQ(filenames.size(),3);
  ASSERT_EQ(filenames[0],jsonFile.getFileDir()+std::string("f0.txt"));
  ASSERT_EQ(filenames[1],jsonFile.getFileDir()+std::string("f1.txt"));
  ASSERT_EQ(filenames[2],jsonFile.getFileDir()+std::string("f2.txt"));
}

BOOST_FIXTURE_TEST_CASE(readArrayArray, JsonFilePaserTestInit){
  std::vector<std::vector<int> > array_array;
  TEST_ASSERT(jsonFile.read("array_array",array_array));
  ASSERT_EQ(array_array.size(), 3);

  ASSERT_EQ(array_array[0].size(), 1);
  ASSERT_EQ(array_array[0][0], 1);

  ASSERT_EQ(array_array[1].size(), 2);
  ASSERT_EQ(array_array[1][0], 2);
  ASSERT_EQ(array_array[1][1], 3);

  ASSERT_EQ(array_array[2].size(), 3);
  ASSERT_EQ(array_array[2][0], 4);
  ASSERT_EQ(array_array[2][1], 5);
  ASSERT_EQ(array_array[2][2], 6);
}

BOOST_AUTO_TEST_SUITE_END()
