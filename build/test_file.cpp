#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include<iostream>
#include <fstream>
using namespace std;
int main(int argv, char* argc[])
{

	srand (time(NULL));
	int elem = atoi(argc[1]);
	ofstream myfile;
	myfile.open ("file1");
	for(int i = 0; i < elem; i++)
	{
		int num = rand() % 200;
		//std::cout << rand() % 200 << std::endl;
	  	myfile << num <<"\n";
	}
	myfile.close();

	myfile.open ("file2");
	for(int i = 0; i < elem; i++)
	{
		int num = rand() % 200;
		//std::cout << rand() % 200 << std::endl;
	  	
	  	myfile << num <<"\n";
	}
	myfile.close();

	myfile.open ("file3");
	for(int i = 0; i < elem; i++)
	{
		int num = rand() % 200;
		//std::cout << rand() % 200 << std::endl;
	  	
	  	myfile << num <<"\n";
	}
	myfile.close();
}
