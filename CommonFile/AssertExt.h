#ifndef ASSERTEXT_H
#define ASSERTEXT_H
/**
 * this is some extetions for the assert command.
 */

#include "Config.h"
#include <string>
#include <stdexcept>
#include <iostream>
#include <boost/lexical_cast.hpp>

#ifndef CUSTOM_DEBUG

# define ASSERTEQ(value_a,value_b) do{}while(0);
# define ASSERTNE(value_a,value_b) do{}while(0);
# define ASSERTGE(value_a,value_b) do{}while(0);
# define ASSERTGT(value_a,value_b) do{}while(0);
# define ASSERTLE(value_a,value_b) do{}while(0);
# define ASSERTLT(value_a,value_b) do{}while(0);
# define ASSERTIN(value_a,min,max) do{}while(0);

#else /* Not NDEBUG.  */

// define a new assert micro to throw expections instead of abort
#define NEWASSERT(e) ((void) ((e) ? 0 : ASSERTTHROW (#e, __FILE__, __LINE__)))
#define ASSERTTHROW( e, file, line ) ( throw std::runtime_error(std::string(file)+":"+boost::lexical_cast<std::string>(line)+": failed assertion "+e))
#ifndef __STRING
#define __STRING(expr) #expr
#endif

// assert "value_a == value_b", if not, print the values and abort.
#define ASSERTEQ(value_a, value_b)										\
  if (value_a != value_b)												\
	{																	\
	  std::cout << __STRING(value_a) <<" = " << value_a<< std::endl;	\
	  std::cout << __STRING(value_b) <<" = " << value_b<< std::endl;	\
      NEWASSERT(value_a == value_b);									\
	}

// assert "value_a != value_b", if not, print the values and abort.
#define ASSERTNE(value_a, value_b)										\
  if (value_a == value_b)												\
	{																	\
	  std::cout << __STRING(value_a) <<" = " << value_a<< std::endl;	\
	  std::cout << __STRING(value_b) <<" = " << value_b<< std::endl;	\
      NEWASSERT(value_a != value_b);									\
	}

// assert "value_a >= value_b", if not, print the values and abort.
#define ASSERTGE(value_a, value_b)										\
  if (value_a < value_b)												\
	{																	\
	  std::cout << __STRING(value_a) <<" = " << value_a<< std::endl;	\
	  std::cout << __STRING(value_b) <<" = " << value_b<< std::endl;	\
      NEWASSERT(value_a >= value_b);									\
	}

// assert "value_a > value_b", if not, print the values and abort.
#define ASSERTGT(value_a, value_b)										\
  if (value_a <= value_b)												\
	{																	\
	  std::cout << __STRING(value_a) <<" = " << value_a<< std::endl;	\
	  std::cout << __STRING(value_b) <<" = " << value_b<< std::endl;	\
      NEWASSERT(value_a > value_b);                                     \
	}

// assert "value_a >= value_b", if not, print the values and abort.
#define ASSERTLE(value_a, value_b)										\
  if (value_a > value_b)												\
	{																	\
	  std::cout << __STRING(value_a) <<" = " << value_a<< std::endl;	\
	  std::cout << __STRING(value_b) <<" = " << value_b<< std::endl;	\
      NEWASSERT(value_a <= value_b);									\
	}

// assert "value_a > value_b", if not, print the values and abort.
#define ASSERTLT(value_a, value_b)										\
  if (value_a >= value_b)												\
	{																	\
	  std::cout << __STRING(value_a) <<" = " << value_a<< std::endl;	\
	  std::cout << __STRING(value_b) <<" = " << value_b<< std::endl;	\
      NEWASSERT(value_a < value_b);                                     \
	}

// assert "value_a in [min,max]", if not, print the values and abort.
#define ASSERTIN(value_a, min, max)                                     \
if ( !(value_a >= min && value_a <= max) )								\
  {																		\
	std::cout << __STRING(value_a) <<" = " << value_a<< std::endl;		\
    std::cout << __STRING(min) <<" = " << min << std::endl;				\
	std::cout << __STRING(max) <<" = " << max << std::endl;				\
    NEWASSERT(value_a >= min && value_a <= max);						\
  }

#endif /* NDEBUG.  */
#endif /* _ASSERTEXT_H_ */
