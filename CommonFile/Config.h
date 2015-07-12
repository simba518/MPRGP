#ifndef CONFIG_H
#define CONFIG_H

#define NAMESPACE COMMON
#define PRJ_BEGIN namespace NAMESPACE {
#define PRJ_END	}
#define USE_PRJ_NAMESPACE using namespace NAMESPACE;

//#define _DEBUG
#ifdef _DEBUG
#define CUSTOM_DEBUG
#endif
#define PROFILE
//#define DOUBLE_PRECISION

#ifdef _MSC_VER
#define ALIGN_16 __declspec(align(16))
#define FORCE_INLINE __forceinline
#else
#define ALIGN_16 __attribute__((aligned (16)))
#define FORCE_INLINE inline
#endif

#include <stdio.h>
//#define NDEBUG
#ifndef NDEBUG
#include <assert.h>
#define ASSERT(x) do{assert((x));}while(0);
#define ASSERT_MSG(x,msg) do{if(!(x)){printf("[ERROR] %s \n",msg);assert(false);}}while(0);

#ifdef _MSC_VER
#define ASSERT_MSGV(x,fmt,...) do{if(!(x)){printf("[ERROR] " fmt " \n",__VA_ARGS__);assert(false);}}while(0);
#else 
#define ASSERT_MSGV(x,fmt,...)
#endif

#else
#ifdef _MSC_VER
#pragma warning(disable:4552)
#pragma warning(disable:4553)
#endif
#define ASSERT(x) do{if(!(x)){exit(EXIT_FAILURE);}}while(0);
#define ASSERT_MSG(x,msg) do{if(!(x)){WARNING(msg);exit(EXIT_FAILURE);}}while(0);
#define ASSERT_MSGV(x,fmt,...) do{if(!(x)){WARNINGV(fmt,__VA_ARGS__);exit(EXIT_FAILURE);}}while(0);
#endif
#define WARNING(msg) do{printf("[WARNING] %s \n",msg);}while(0);

#ifdef _MSC_VER
#define WARNINGV(fmt,...) do{printf("[WARNING] " fmt " \n",__VA_ARGS__);}while(0);
#else
#define WARNINGV(fmt,...) 
#endif

#define INFO(msg) do{printf("[INFO] %s \n",msg);}while(0);
#define INFOV(fmt,...) do{printf("[INFO] " fmt " \n",__VA_ARGS__);}while(0);
#define NOTIFY_MSG(msg) do{printf("[NOTIFY] %s \n",msg);}while(0); getchar();
#define NOTIFY_MSGV(fmt,...) do{printf("[NOTIFY] " fmt " \n",__VA_ARGS__);}while(0); getchar();

//OpenMP only support signed variable as index
#include <stdint.h>
#if defined(WIN32) || defined(__CUDACC__)
typedef int64_t sizeType;
#else
//#define sscanf_s sscanf
//#define sprintf_s sprintf
//typedef long long sizeType;
typedef int64_t sizeType;
#endif

typedef int vtkSizeType;
#define INVALID sizeType(-1)

#endif
