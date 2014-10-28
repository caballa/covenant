project(covenant)

cmake_minimum_required (VERSION 2.8.11)

if (CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message (FATAL_ERROR "In-source builds are not allowed. Please clean your source tree and try again.")
endif()

# Add path for custom modules
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

#set(CMAKE_CXX_FLAGS "-g -O3 -Wall -Wno-non-virtual-dtor -Wno-unused-function -Wno-deprecated -ffast-math")
set(CMAKE_CXX_FLAGS "-g -O3 -Wall")

set (CUSTOM_BOOST_ROOT "" CACHE PATH "Path to custom boost installation.")
if (CUSTOM_BOOST_ROOT)
  set (BOOST_ROOT ${CUSTOM_BOOST_ROOT})
  set (Boost_NO_SYSTEM_PATHS "ON")
endif()

set (Boost_USE_STATIC_LIBS ON)
find_package (Boost 1.55.0 COMPONENTS program_options REQUIRED)
IF (Boost_FOUND)
  include_directories (${Boost_INCLUDE_DIRS})
endif ()

#add_definitions(-std=c++11)
#add_definitions(-Wno-redeclared-class-member -Wno-sometimes-uninitialized)
#add_definitions(-Wno-unused-local-typedefs)

## This should be enabled until stable version
add_definitions(-DSANITY_CHECKS)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
add_subdirectory(lib)
add_subdirectory(tools)





