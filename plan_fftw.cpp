#include "header.h"
#include <blitz/array.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

using namespace blitz;
using namespace std;

#define FFTW_PLAN_TYPE FFTW_MEASURE

void GB_system::plan_fftw() 
{
 	//nthreads = 4;
   	//cout << "FFT plan initialized using " << nthreads << " threads." << endl;
 	//fftw_plan_with_nthreads(nthreads);
	
	fftw_p1 = fftw_plan_dft_3d(Nz, Ny, Nx, reinterpret_cast<fftw_complex*>(psi1.data()) , reinterpret_cast<fftw_complex*>(F1.data()) , FFTW_FORWARD, FFTW_PLAN_TYPE);
	fftw_p2 = fftw_plan_dft_3d(Nz, Ny, Nx, reinterpret_cast<fftw_complex*>(psi2.data()) , reinterpret_cast<fftw_complex*>(F2.data()) , FFTW_FORWARD, FFTW_PLAN_TYPE);
	fftw_p3 = fftw_plan_dft_3d(Nz, Ny, Nx, reinterpret_cast<fftw_complex*>(psi3.data()) , reinterpret_cast<fftw_complex*>(F3.data()) , FFTW_FORWARD, FFTW_PLAN_TYPE);
	fftw_p4 = fftw_plan_dft_3d(Nz, Ny, Nx, reinterpret_cast<fftw_complex*>(psi4.data()) , reinterpret_cast<fftw_complex*>(F4.data()) , FFTW_FORWARD, FFTW_PLAN_TYPE);
	ifftw_F1 = fftw_plan_dft_3d(Nz, Ny, Nx, reinterpret_cast<fftw_complex*>(F1.data()) , reinterpret_cast<fftw_complex*>(psi1.data()) , FFTW_BACKWARD, FFTW_PLAN_TYPE);
	ifftw_F_mu = fftw_plan_dft_3d(Nz, Ny, Nx, reinterpret_cast<fftw_complex*>(F_mu.data()) , reinterpret_cast<fftw_complex*>(mu.data()) , FFTW_BACKWARD, FFTW_PLAN_TYPE);
	ifftw_FE_term = fftw_plan_dft_3d(Nz, Ny, Nx, reinterpret_cast<fftw_complex*>(FE_term.data()) , reinterpret_cast<fftw_complex*>(E_term.data()) , FFTW_BACKWARD, FFTW_PLAN_TYPE);
}
