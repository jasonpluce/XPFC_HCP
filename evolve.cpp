#include "header.h"
#include <blitz/array.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <fstream>
#include <time.h>
#include <string>

using namespace blitz;
using namespace std;

void GB_system::evolve() 
{	
	Array<complex<double>,3> term1(Nz,Ny,Nx), term2(Nz,Ny,Nx), p2(Nz,Ny,Nx), p3(Nz,Ny,Nx), p_temp(Nz,Ny,Nx), filtered_p(Nz,Ny,Nx);
	Array<complex<double>,3> E_ext(Nz,Ny,Nx), E_int(Nz,Ny,Nx), LED(Nz,Ny,Nx), F1_temp(Nz,Ny,Nx), fft_mu_0(Nz,Ny,Nx);
 	Range F2_dealias_Nz(2*Nz/3,Nz-1), F2_dealias_Ny(2*Ny/3,Ny-1), F2_dealias_Nx(2*Nx/3,Nx-1), F3_dealias_Nz(Nz/2,Nz-1), F3_dealias_Ny(Ny/2,Ny-1), F3_dealias_Nx(Nx/2,Nx-1);

	int count = 1, E_step = 100, n = 0;
 	double dE = 1e15, E_old = 0.0, Energy = 0.0, F1_0 = 0.0, F_mu_0 = 0.0, psi_mean_from_F1_0 = 0.0, mu_mean_from_F_mu_0 = 0.0;
	double psi_mean_initial = 0.0, psi_mean = 0.0, mu_mean = 0.0;
	double mu_diff = 0.0, conv_thresh = 5e-8, N_tot = Nx*Ny*Nz;
 	complex<double> conserved_density_val = 0.0, fft_mu_0_val = 0.0;

	cout << setprecision(14);
	n = 1000000;
	fft_mu_0 = 0.0;
	
	//write_out(real(psi1),"psi" + to_string((int)target_angle) + "_0");//output vtk file for IC
 	//write_out(real(psi1),"psi" + to_string((int)s_val) + "_0"); // output vtk file for strain testing	
	
	// UC tests for fft_mu_0, conv_thresh = 5e-8
		// dx = 0.1, fft_mu_0 = -0.066710385076912
		// dx = 0.09, fft_mu_0 = -0.066721871454396
		// dx = 0.08, fft_mu_0 = -0.066705626761045
	
	fft_mu_0_val = -0.066721871454252*N_tot;
	fft_mu_0(0,0,0) = fft_mu_0_val;

	psi_mean_initial = mean(real(psi1));
	fftw_execute(fftw_p1);
	
	conserved_density_val = n_bar*N_tot;	
	psi2 = psi1*psi1;
 	psi3 = psi1*psi2;

	cout << "Evolving system" << endl;
	clock_t t;
	t = clock();

 	for (int iter=0; iter<n; iter++)
 	{

 		fftw_execute(fftw_p2);
 		fftw_execute(fftw_p3);

		//F2(F2_dealias_Nx,F2_dealias_Ny,F2_dealias_Nz) = 0.0; // dealias psi^2 by setting last 1/3 of wavenumber values to 0; CHECK COMPLEX TERMS
		//F3(F3_dealias_Nx,F3_dealias_Ny,F3_dealias_Nz) = 0.0; // dealias psi^2 by setting last 1/2 of wavenumber values to 0;

 		// term1 = F1 + (F2/2.0 - F3/3.0)*ks*dt; // locally conserved dynamics
 		// term2 = 1.0 + (1.0 - C)*ks*dt;
 		
		term1 = F1 + (F2/2.0 - F3/3.0 + fft_mu_0)*dt; //non-conserved dynamics with constant mu
        term2 = 1.0 + (1.0 - C)*dt; 
		F1 = term1/term2;
		//F1(0,0,0) = conserved_density_val; // used to conserve mass
		
		fftw_execute(ifftw_F1);

 		psi1 = psi1/N_tot;// rescale psi
		psi2 = psi1*psi1;
 		psi3 = psi1*psi2;

 		if (iter%E_step == 0)
 		{
 			FE_term = C*F1;
 			fftw_execute(ifftw_FE_term);
 			E_term = E_term/N_tot;
 			
 			E_int = psi2/2.0 - psi3/6.0 + (psi3*psi1)/12.0;
 			E_ext = -0.5*psi1*E_term;
        	
			F_mu = F1 - F2/2.0 + F3/3.0 - FE_term; // FFT of chemical potential, e.g. FFT(mu(r)),
			cout << "Fixed mu value: "  <<	fft_mu_0_val << endl; // value set for constant chemical potential
			cout << "F_mu_0: " << F_mu(0,0,0) << endl << endl;// value of F_mu zero mode, should converge to fft_mu_0_val when mass is not being conserved
			fftw_execute(ifftw_F_mu); // returns mu 
        	mu = mu/N_tot;
			// mu = psi1 - psi2/2.0 + psi3/3.0 - E_term; // alternative definiton of mu - should be identical to ifft(F_mu) 

			mu_diff = max(real(mu)) - min(real(mu));// determines largest gradient between any two points in the chemical potential - used for convergence
			LED = real(E_int + E_ext); // Local energy density landscape
        	Energy = sum(real(E_int + E_ext))*dx*dy*dz;

			F1_0 = real(F1(0,0,0));
			F_mu_0 = real(F_mu(0,0,0));
			
			psi_mean  = mean(real(psi1));
        	mu_mean = mean(real(mu));
			
			psi_mean_from_F1_0  = F1_0/N_tot;
			mu_mean_from_F_mu_0 = F_mu_0/N_tot;
			
			ED = Energy/(Lx*Ly*Lz);
        	dE = (Energy-E_old)/E_step;
        	E_old = Energy;

        	cout << "Time step: " << iter << endl;
			cout << "dE: " << dE << endl;
        	cout << "mu diff: " << mu_diff << endl;
			//cout << "psi_mean_from_F1_0: " << psi_mean_from_F1_0 << " psi_mean: " << psi_mean  << " mu_mean_from_F_mu_0: " << mu_mean_from_F_mu_0 << " mu_mean: " << mu_mean<< endl << endl;
			cout << "psi_mean_from_F1_0: " << psi_mean_from_F1_0 << " mu_mean_from_F_mu_0: " << mu_mean_from_F_mu_0 << endl << endl;
		
        	if (std::isnan(Energy))
        	{
        		cout << "Solution is unstable - evolution terminating." << endl;
        		//cout << "Time step: " << i << " Energy: " << Energy << " Energy density: " << Energy/(Lx*Ly*Lz) << " dE: " << dE << endl;
        		break;
        		//exit (EXIT_FAILURE);
        	}
        
        	if (mu_diff < conv_thresh)
        	{
        		cout << "Solution has converged - evolution terminating." << endl;
        		//cout << "Time step: " << i << " Energy: " << Energy << " Energy density: " << Energy/(Lx*Ly*Lz) << " dE: " << dE << endl;
        		break;
        		//exit (EXIT_SUCCESS);
			}

			// if (abs(dE) < conv_thresh)
        	// {
        	// 	cout << "Solution has converged - evolution terminating." << endl;
        	// 	//cout << "Time step: " << i << " Energy: " << Energy << " Energy density: " << Energy/(Lx*Ly*Lz) << " dE: " << dE << endl;
        	// 	break;
        	// 	//exit (EXIT_SUCCESS);
        	// }

        	else
        	{
        		//cout << "i: " << i << endl;
        		//write_out(real(psi1),"psi" + to_string((int)target_angle) + "_" + to_string(count));
        		//write_out(real(psi1),"psi" + to_string((int)s_val) + "_" + to_string(count)); 
        	}
        	count++;
 		} 	
	}
	// cout << "Solution has timed out - evolution terminating." << endl;
	cout << "Lx: " << Lx << '\t' << "Ly: " << Ly << '\t' << "Lz: " << Lz <<endl;
	cout << "dx: " << dx << '\t' << "dy: " << dy << '\t' << "dz: " << dz <<endl;
	cout << "angle: " << theta*RAD2DEG << endl;
	cout << "ED: " << ED << endl;
	cout << "dE: " << dE << endl;
	cout << "mu diff: " << mu_diff << endl;
	//cout << "psi_mean_from_F1_0: " << psi_mean_from_F1_0  << " psi_mean: " << psi_mean  << " mu_mean_from_F_mu_0: " << mu_mean_from_F_mu_0 << " mu_mean: " << mu_mean << endl << endl;
	cout << "psi_mean_from_F1_0 : " << psi_mean_from_F1_0 << " mu_mean_from_F_mu_0: " << mu_mean_from_F_mu_0 << endl << endl;

	//cout << "Energy: " << Energy << " Energy density: " << ED << " dE: " << dE/E_step << " fft_mu_0: " << fft_mu_0_val << endl;
	//write_out(real(psi1),"psi" + to_string((int)target_angle));
	//write_out(real(psi1),"psi" + to_string((int)s_val));
	//write_out(real(F1),"F1");
	//write_out(real(FE_term),"FE_term");
	//cout << "F1_mean: " << mean(real(F1)) << endl;
	cout << endl;

	// Output results to .txt file - BAD!!! update to better format then delete this section
	// ofstream GBE_data;
	// GBE_data.open("GBE_data.txt");
	// GBE_data.open("GBE_data.txt", ofstream::app);
	// GBE_data.precision(16);
	// GBE_data.setf(std::ios::fixed, std:: ios::floatfield); // floatfield set to fixed
  	// GBE_data << "angle: " << RAD2DEG*theta << endl;
  	// GBE_data << "Energy: " << Energy << endl;
  	// GBE_data << "Energy Density: "  << Energy/(Lx*Ly*Lz) << endl;
  	// GBE_data << "Lx: " << Lx << endl;
  	// GBE_data << "Ly: " << Ly << endl;
  	// GBE_data << "Lz: " << Lz << endl;
  	// GBE_data << "Nx: " << Nx << endl;
  	// GBE_data << "Ny: " << Ny << endl;
  	// GBE_data << "Nz: " << Nz << endl<< endl;
	// GBE_data.close();
	
	t = clock() - t;
	
	//if (s_val <= 20) {GB_system system(target_angle,s_val);}
		
	
	//cout << endl << "Energy density: " << Energy/(Lx*Ly*Lz) << endl;
	cout << "Time elapsed: " << (double)t/CLOCKS_PER_SEC << " seconds" << endl << endl;

//	write_out(real(psi1),"psi2.vtk"); 
// 	write_out(psi2,"psi_squared.txt");
// 	write_out(psi3,"psi_cubed.txt");
// 	write_out(term1,"term1.txt");
// 	write_out(term2,"term2.txt");
// 	write_out(psi1_new_hat,"F1.txt");
	
}
