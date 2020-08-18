#include <iostream>
#include <string>

#include "header.h"

// Include Blitz++ (This makes it nice to deal with multi-dimensional arrays)
#include <blitz/array.h>
#include <string>
//#include <omp.h>

using namespace blitz;
using namespace std;

//
// constructor for unit_cell class
//

XPFC_system::XPFC_system(string input_file) 
{
	
	//define parameters (values will eventually be assigned via input file)
	
  	//dx_base = 0.05;
  	//dy_base = 0.05;
  	//dz_base = 0.05;
	
 	dx_base = 0.1;
	dy_base = 0.1;
 	dz_base = 0.1;
	
	dt = 0.10;
	
	dim = 3; 
	a = 1.0;
	c = sqrt(8.0/3.0)*a;
	q = 2.0*PI/(a); 
	
	k1 = 4.0*PI/(a*sqrt(3.0)); 
	k2 = 2.0*PI*sqrt(41.0/24.0)/a;
	
	alpha1 = 1.0;
	alpha2 = alpha1;
	alpha3 = alpha1;
	
	sigma = 0.00;
	psi_mean = 0.0; 
	n_bar = 0.0;
	//n_bar = 0.0025;
	//n_bar = 0.0319;
	mass = 0.0;
 	
	beta1 = 8.0; 
	beta2 = 6.0;
	beta3 = 4.0;
	
	rho1 = 2.0; 
	rho2 = 2.0;
	rho3 = 2.0;
	
	//2D parameters
	sigma = 0.8/(2.0*PI); 
	n_bar = 0.14;
	rho1 = 1.0;
	beta1 = 4.0;
	
	chi = 0.5;
	nu = 1.0/6.0;
	epsi = 1.0/12.0;
	
// 	target_angle = angle;
// 	s_val = size_val;
	//target_angle = 20.0;
	da_thresh = 0.5;
	//IC_spacing = 0;
	IC_spacing = 5;
	
	corr_type = "HEXAGONAL";
	IC_type = "HEXAGONAL";
	
	//corr_type = "BCC";
	//IC_type = "BCC_UC";
	
	//corr_type = "HCP";
	//IC_type = "HCP_UC";
	//IC_type = "HCP_YROT";
	//IC_type = "HCP_ZROT";
	//IC_type = "FCC_UC";
	//IC_type = "Liquid";
	
//	nthreads = omp_get_max_threads();

	get_size(); // determine values for strain free system size;
	
// resize initialized arrays
	psi1.resize(Nz,Ny,Nx); psi2.resize(Nz,Ny,Nx); psi3.resize(Nz,Ny,Nx); psi4.resize(Nz,Ny,Nx); 
	F1.resize(Nz,Ny,Nx); F2.resize(Nz,Ny,Nx); F3.resize(Nz,Ny,Nx); F4.resize(Nz,Ny,Nx);
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
