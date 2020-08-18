#ifndef HEADER_H
#define HEADER_H

#include <math.h>
#include <complex>

#include "fftw3.h"
#include <blitz/array.h>
#include <string>

using namespace blitz;
using namespace std;

// define constant global variables
const double PI = acos(-1); // define PI
const double DEG2RAD = (PI/180.0); // used to convert from degrees to radians
const double RAD2DEG = (180.0/PI); // used to convert from radians to degrees

class GB_system
{	
	void get_size(); // determine strain free system size 
	void plan_fftw(); // make plans for FFTW
	void make_C();
	
	public:
		
		//GB_system(double angle, double size_val);
		GB_system(string input_file);
		
		// declare public variables
		
		int nthreads; // number of threads used by FFTW
		
		double dx_base, dy_base, dz_base;// Guide for grid spacing values. Actually grid spacing may be slightly larger or smaller depending on strain free system size
		double dt; // time step
		
		int dim; // dimensionality of the system 2 = 2D, 3 = 3D (1D not supported)

		// structural parameters
		double a,c; // unit cell lattice parameters
		double q; // inter-planar spacing for IC
		double k1,k2,k3; // location of the peak(s) of the correlation function in Fourier space
		double alpha1,alpha2,alpha3; // width of the peak(s) of the correlation function
		double mass; // mass of the system;
		double psi_mean; // average value of the initialized order parameter
		double n_bar; // desired average value for the order parameter (i.e. average normalized density of the system)
		double ED; // energy density
		double chi, nu, epsi; //fitting parameters for the XPFC free energy 
	
		double sigma; // weighting term of XPFC Debye-waller factor (determines height of all peaks of the correlation function - goes like 1/exp(k)
		double beta1,beta2,beta3; // planar symmetry - phenomenological
		double rho1,rho2,rho3; // planar density - phenomenological
	
		// integer value for determining unstrained system size
		int l; // determines Lx and Ly (Lz to be added)

		// grain mis-orientation (expressed as a symmetric tilt angle) and angular resolution threshold(i.e. margin of acceptable error in target_angle)
		double target_angle;// target symmetric tilt angle between the grains, in degrees
		double theta;// strain free angle between grains, in radians (difference between target_angle and RAD2DEG*theta is bounded by da_thresh)
		double da_thresh; // tolerance for thet difference between target_angle and RAD2DEG*theta, in degrees
						  // smaller values for 'da_thresh' generally result in a larger system size for a given target_angle)
		
		//double angle_start, angle_stop; // start and and stop angle for grain boundary anisotropy simulations 
		
		// define width of liquid region surrounding the grains for the IC
		double IC_spacing;

		//define variables
		double Lx, Ly, Lz; //System size in each direction
		double Nx, Ny, Nz; //Number of grid points in each directions
		double dx,dy, dz; // grid spacing in each direction
		
		double s_val; // baseline value used to determine system size
		
		string corr_type, IC_type; //flag that determines the correlation function, currently supports "HEXAGONAL" and "HCP"
		
		Array<complex<double>,3> psi1; //order parameter used to simulate atomic density
		Array<complex<double>,3> psi2, psi3, psi4; // psi^2,psi^3, and psi^2 + psi^3
		Array<complex<double>,3> mu; // chemical potential
		Array<complex<double>,3> F_mu; // chemical potential
		Array<complex<double>,3> F1; // Fourier transform of psi1
		Array<complex<double>,3> F2, F3, F4; // Fourier transform of psi^2 and psi^3
		Array<complex<double>,3> E_term, FE_term; // used to calculate the energy of the system
		
		Array<complex<double>,3> ks; // k^2 - used in XPFC dynamics 
		Array<complex<double>,3> C; // XPFC correlation function 

		fftw_plan fftw_p1, fftw_p2, fftw_p3, fftw_p4, ifftw_F1, ifftw_FE_term, ifftw_F_mu; // plans for use with FFTW
	
		GB_system(); //constructor
		void set_IC(); // build initial conditions using OMA
		void set_angle(); // set symmetric tilt angle of grains
		void evolve(); // evolve the system to equilibrium conditions
		void write_out(Array<double,3> array, string filename); // write output to file	
		void write_out(Array<complex<double>,3> array, string filename); // write output to file
		void write_out_smooth(Array<double,3> array, string filename); // write output to file	
		void write_out_smooth(Array<complex<double>,3> array, string filename); // write output to file
		void test_strain();	
};

#endif
	
