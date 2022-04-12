cmake_minimum_required (VERSION 3.11)

project(LigeroLink CXX)

set(Boost_USE_STATIC_LIBS     ON)
set(Boost_USE_MULTITHREADED   ON)

find_package(Boost
  COMPONENTS program_options system log exception serialization
  REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include ${Boost_INCLUDE_DIR})

add_compile_options(-Wall -std=c++20 -fvisibility=hidden -stdlib=libc++)

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  add_compile_options(-O3 -flto=thin)
  add_link_options(-flto=thin)
endif()

# ------------------------------------------------------------

add_library(portable_binary_archive
  STATIC
  include/ligero/util/boost/portable_binary_iarchive.cpp
  include/ligero/util/boost/portable_binary_oarchive.cpp)
target_link_libraries(portable_binary_archive ${Boost_LIBRARIES} c++ c++abi)
