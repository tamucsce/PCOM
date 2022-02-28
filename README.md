Installation

GMP
Install gmp-6.2.1 follow official instruction or:
  cd gmp
  ./configure
  make
  make check
  make install

Boost
Install boost_1_76_0
 

Build Instruction
To run the PSI protocal between 2 parties:

open 2 terminals

both:

 cd /MultipartyPSI/build 
compile:

 make
run:

 terminal 1: ./src/mpsi -s -f file1

 terminal 2: ./src/mpsi -c -f file2
