cmake_minimum_required (VERSION 2.8.11)
project (multiparty-psi)
set(NAME "multiparty-psi")



#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)

#add_compile_definitions(BN_SUPPORT_SNARK MIE_ATE_USE_GMP)






set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} cmake/)
include_directories(${CMAKE_SOURCE_DIR})
include_directories("../depends/ate-pairing/include/")
include_directories("../depends/ate-pairing/tests/")
link_directories("../depends/ate-pairing/lib/")


include_directories("../depends/KZG_ext-main/include/")


set (CMAKE_C_FLAGS "-pthread  -march=native -O3 -maes")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg -ggdb")
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "-DBN_SUPPORT_SNARK")

add_executable(
    mpsi 
    mpsi.cpp
)
target_link_libraries(
    mpsi

    gmp
    gmpxx
    zm

)

