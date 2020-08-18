//========================================================================================
// main.cpp
//----------------------------------------------------------------------------------------


//#include <iostream>
//#include <stdio.h>
#include <string>
//#include <sstream>
//#include <math.h>
#include "fftw3.h"
#include "header.h"
#include <blitz/array.h>
#include <complex.h>
#include <fstream>
//#include <omp.h>

using namespace blitz;
using namespace std;

// Program starts here
int main ()
{	
	string input_file;
	input_file = "test.dat";
 	int angle_start, angle_stop,sv_start,sv_stop,di,dv;
 	
 	Array<double,1> angle_vals(34);
 	//angle_vals = 0.1, 5.9, 8.9, 11.7, 13.1, 17.2, 21.1, 22.4, 24.9, 27.3, 28.4, 31.7, 34.8, 35.8, 37.7, 39.5, 40.3, 42.8, 45.1, 45.9, 47.2, 48.5, 51.0, 54.3, 57.1, 59.5, 61.6, 63.5, 65.2, 68.0, 70.2, 74.9, 78.6, 90.0; 
 	angle_start = 0;
 	angle_stop = 0;
 	sv_start = 0;
 	sv_stop = 0;
 	di = 1;
 	dv = 1;
	int fftw_init_threads(void); // initialize FFTW multithreading
	
 	for (int i = angle_start; i<=angle_stop; i=i+di)
 	//for (int i = sv_start; i<=sv_stop; i=i+dv)
 	{
		//GB_system system(angle_vals(i),0.0);
		GB_system system(input_file);
  		system.set_IC();
  		system.evolve();
  		//system.test_strain();
 	}
 	
 	//void fftw_cleanup_threads(void); // clean up threads
 	
// string input_file;
//	const char * input_file;
//	input_file="input.dat";
	
// 	cout << endl << "Starting FCC section" << endl;
// 	unit_cell fcc_object(input_file,"FCC");
// 	fcc_object.Loop_Params();
	
//	cout << endl << "Starting BCC section" << endl;
//	unit_cell bcc_object(input_file,"BCC");
//	bcc_object.Loop_Params();
	
// 	cout << endl << "Starting diamond section" << endl;
// 	unit_cell diamond_object(input_file,"diamond");
// 	diamond_object.Loop_Params(); 
 	
//	cout << endl << "Program complete" << endl;
// 	bcc_object.Write_To_Log();


} //end of main
