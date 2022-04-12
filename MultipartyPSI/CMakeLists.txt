cmake_minimum_required(VERSION 3.5)

PROJECT (MultipartyPSI CXX)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(Boost_USE_STATIC_LIBS     ON)
set(Boost_USE_MULTITHREADED   ON)

find_package(GMP REQUIRED)
find_package(Boost
  COMPONENTS program_options system log exception serialization
  REQUIRED)

# Include headers
# ------------------------------------------------------------
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/depends/xbyak/xbyak/)
include_directories(${PROJECT_SOURCE_DIR}/depends/ate-pairing/include/)
include_directories(${PROJECT_SOURCE_DIR}/depends/ate-pairing/tests/)
include_directories(${PROJECT_SOURCE_DIR}/depends/KZG_ext-main/include/)
# LigeroLink
include_directories(${PROJECT_SOURCE_DIR}/depends/LigeroLink/include)
# Third-party
include_directories(${GMP_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIR})

# Compile options
# ------------------------------------------------------------
if(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  add_compile_options(-O3 -march=native)
endif()

add_compile_options(-fvisibility=hidden)

# Build libraries
# ------------------------------------------------------------

# Build LigeroLink example
add_subdirectory(${PROJECT_SOURCE_DIR}/depends/LigeroLink/  CMAKE_CURRENT_BINARY_DIR)
add_executable(demo_add examples/Add.cpp)
target_compile_options(demo_add PUBLIC -std=c++20 -stdlib=libc++)
target_link_libraries(demo_add ${Boost_LIBRARIES} portable_binary_archive c++ c++abi)

# Build MultipartyPSI/src
add_compile_definitions(BN_SUPPORT_SNARK MIE_ATE_USE_GMP)
add_library(zm STATIC
  ${PROJECT_SOURCE_DIR}/depends/ate-pairing/src/zm.cpp
  ${PROJECT_SOURCE_DIR}/depends/ate-pairing/src/zm2.cpp)

add_subdirectory(src)
set(CMAKE_CXX_FLAGS "-DBN_SUPPORT_SNARK")
