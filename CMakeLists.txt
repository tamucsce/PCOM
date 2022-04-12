cmake_minimum_required (VERSION 3.5)
project (multiparty-psi)
set(NAME "multiparty-psi")


if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  add_compile_options(-std=c++11 -pg -ggdb)
endif()

if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
  add_compile_options(-O3 -march=native -maes)
endif()

#add_executable(
#    mpsi 
#    mpsi.cpp
#)
#target_link_libraries(mpsi GMP::gmp GMP::gmpxx zm)

add_executable(semihonest semihonest.cpp)
target_compile_options(semihonest PUBLIC -std=c++20 -fvisibility=hidden -stdlib=libc++)
target_link_libraries(semihonest GMP::gmp GMP::gmpxx portable_binary_archive ${Boost_LIBRARIES} zm c++ c++abi)
