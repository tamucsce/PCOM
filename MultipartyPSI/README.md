# MultipartyPSI

## Building the Examples

In order to build the library and examples, you'll need `Boost >= 1.74` and a C++20 compiler.


``` bash
mkdir build
cd build/
cmake -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_BUILD_TYPE=Release ../examples
make
```

The `CMAKE_CXX_COMPILER` flag is optional if your system compiler support `-std=c++20`.


## MultipartyPSI

Here everything the current model can do:

1. malicious part of MultipartyPSI (multipoints Eval)
2. polynomial multiplication, division
3. pairing
4. roots -> poly coeffs (interpolation)



in progress (Wenxuan & Yupeng):

1. Use REAL encryption instead of commitment 
2. Add Fiat-Shamir 
3. more optimizing




Usage:

$ cd src

$ make

$ ./mpsi \<log degree of polynomial\>


