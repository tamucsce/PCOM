Here is the ami image for our AWS instance:
N.Virginia

ami-0fa60e2ac8161c6b0

It is already set up. 

To run the PSI protocal between 2 parties:

1) open 2 terminals

    both:
	cd /MultipartyPSI/build 

    compile:
	make

    run:
	terminal1: ./src/mpsi -s -f file1
	terminal2: ./src/mpsi -c -f file2






*NOTE: -s denotes the P1 party, -c denotes other parties. Each party's set is in each corresponding file. To generate the file, there is a code called "test_file.cpp" that generates it (Or you can import them yourself)
	If you want to test it on 3 or more parties, just open more terminals.
	Due to limited time, this code is implemented to only support polynomials that has size exactly to the power of 2. So for now only play with it with input files with 1 less the power of 2 amount of numbers (EX: 7, 63, 1023).
	Sorry for inconvenience.
	




