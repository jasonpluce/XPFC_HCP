#include <iostream>
#include <string>

#include "header.h"

// Include Blitz++ (This makes it nice to deal with multi-dimensional arrays)
#include <blitz/array.h>
#include <string>
//#include <omp.h>

using namespace blitz;
using namespace std;

//GB_system::GB_system(double angle, double size_val) 
GB_system::GB_system(string input_file) 
{
	
	//define parameters (values will eventually be assigned via input file)
	
	//string input_file;
	//const char * input_file;
	//input_file="HCP_infile.txt";
	
	//2D parameters
	// 	sigma = 0.8/(2.0*PI); 
	// 	n_bar = 0.14;
	// 	rho1 = 1.0;
	// 	beta1 = 4.0;

	cout << endl << "input file: "  << input_file << endl;

	da_thresh = 0.5;
	target_angle = 11.0;

 	dt = 1.0;
	
	dx_base = 0.1;
	dy_base = dx_base;
 	dz_base = dx_base;
	
	dim = 3; 
	
	a = 1.0;
	c = sqrt(8.0/3.0)*a;
	q = 2.0*PI/a;
	
	alpha1 = 1.0;
	alpha2 = alpha1;
	alpha3 = alpha1;
	
	sigma = 0.00;
	n_bar = 0.0;
	mass = 0.0;
 	
	beta1 = 8.0; 
	beta2 = 6.0;
	beta3 = 4.0;
	
	rho1 = 2.0; 
	rho2 = rho1;
	rho3 = rho1;
	
	chi = 0.5;
	nu = 1.0/6.0;
	epsi = 1.0/12.0;
	
	//corr_type = "HEXAGONAL";
	//IC_type = "HEXAGONAL";
	
	//corr_type = "BCC";
	//IC_type = "BCC_UC";
	
	corr_type = "HCP";
	IC_type = "HCP_UC";
	//IC_type = "HCP_YROT";
	//IC_type = "HCP_ZROT";

	//IC_type = "FCC_UC";
	//IC_type = "BCC_UC";
	//IC_type = "Liquid";
	//IC_type = "Columnar";
	//IC_type = "Spherical";
	
//	nthreads = omp_get_max_threads();

	get_size(); // determine values for strain free system size;
	
// resize initialized arrays
	psi1.resize(Nz,Ny,Nx); psi2.resize(Nz,Ny,Nx); psi3.resize(Nz,Ny,Nx); psi4.resize(Nz,Ny,Nx); 
	F1.resize(Nz,Ny,Nx); F2.resize(Nz,Ny,Nx); F3.resize(Nz,Ny,Nx); F4.resize(Nz,Ny,Nx);
	F_mu.resize(Nz,Ny,Nx); mu.resize(Nz,Ny,Nx); 
	E_term.resize(Nz,Ny,Nx); FE_term.resize(Nz,Ny,Nx);
	C.resize(Nz,Ny,Nx); ks.resize(Nz,Ny,Nx);

// make FFTW plans *** IMPORTANT *** - array values are (usually) wiped when FFT plans are made - therefore the FFT plan should be made before initializing array values 
 	plan_fftw();
// make correlation function
 	make_C();
		
	//this->test_strain();
		
} //End unit_cell constructor

// destructor for unit_cell class (cleans up memory)
//unit_cell::~unit_cell () {
	//Remember to delete any dynamic memory allocated
	//delete n;
	//Destroy Fourier Transform (may not want to do this if we use this function every time iteration...)
	//fftw_destroy_plan(fftw_n);
	//fftw_destroy_plan(fftw_n2);
	//fftw_destroy_plan(fftw_n3);
	//fftw_destroy_plan(ifftw_mu_hat);
	//fftw_destroy_plan(ifftw_corr_hat); 
//} // end of unit_cell destructor
