#include "header.h"
#include <blitz/array.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace blitz;
using namespace std;

void GB_system::test_strain() 
{
	cout << endl << "Beginning unit cell strain test" << endl;

	double Lx_base = Lx;
	double Ly_base = Ly;
	double Lz_base = Lz;
	
	double dx_base = dx;
	double dy_base = dy;
	double dz_base = dz;
	
	//double rho_0 = 2.0;
	
	Array<double,1> ED_vals_Lx;
	Array<double,1> ED_vals_Ly;
	Array<double,1> ED_vals_Lz;
	Array<double,1> ED_vals_V;

	//int dm = 70, dp = 95;
	int dm = -10, dp = 10;
	double inc = 0.001;
	int count = 0;

	ED_vals_Lx.resize(dp-dm+1);
	ED_vals_Ly.resize(dp-dm+1);
	ED_vals_Lz.resize(dp-dm+1);
	ED_vals_V.resize(dp-dm+1);
	
	//psi1.resize(Nz,Ny,Nx); psi2.resize(Nz,Ny,Nx); psi3.resize(Nz,Ny,Nx); 
	//F1.resize(Nz,Ny,Nx); F2.resize(Nz,Ny,Nx); F3.resize(Nz,Ny,Nx);
	//E_term.resize(Nz,Ny,Nx); FE_term.resize(Nz,Ny,Nx);
	//C.resize(Nz,Ny,Nx); ks.resize(Nz,Ny,Nx);
	
	//n_bar = 0.000;
	//n_bar = 0.0319;
	//n_bar = 0.031874385467731;
	double n_base = n_bar;
	double test_mass;
	
	//this->set_IC();
	//double base_mass = sum(real(psi1))*dx*dy*dz;
	double base_mass = sum(real(psi1)+1.0)*dx*dy*dz;
	cout << "dx: " << dx << endl;
	cout << "base mass: " << base_mass << endl << endl;

	for (int i=dm; i<=dp; i++)
	{
		Lx = Lx_base + inc*Lx_base*i;
		dx = Lx/(Nx);
	 	n_bar = (n_base*Lx_base + (Lx_base-Lx))/Lx;
	 	
		//plan_fftw();
		make_C();
		//this->set_IC();
		
		cout << "n_bar: " << n_bar << endl;
		//psi1 = psi1+n1;
		//test_mass = sum(real(psi1))*dx*dy*dz;
		test_mass = sum(real(psi1)+1.0)*dx*dy*dz;
 		cout << "test mass: " << test_mass << endl;
 		psi_mean = mean(real(psi1));
 		cout << "psi mean: " << psi_mean << endl;
		this->evolve();
		ED_vals_Lx(count) = ED;
		count++;
		write_out(real(psi1),"psi" + to_string(count));
	}
// 	
// 	Lx = Lx_base;
// 	dx = dx_base;
// 	count = 0;
// 	
// 	for (int i=dm; i<=dp; i++)
// 	{
// 		Ly = Ly_base + inc*Ly_base*i;
// 		dy = Ly/(Ny-1);
// 		n_bar = (n_base*Ly_base + (Ly_base-Ly))/Ly;
// 		
// 		plan_fftw();
// 		make_C();
// 		this->set_IC();
// 		this->evolve();
// 		ED_vals_Ly(count) = ED;
// 		count++;
// 	}
// 	
// 	Ly = Ly_base;
// 	dy = dy_base;
// 	count = 0;
// 	
// 	for (int i=dm; i<=dp; i++)
// 	{	
// 		Lz = Lz_base + inc*Lz_base*i;
// 		dz = Lz/(Nz-1);
// 		n_bar = (n_base*Lz_base + (Lz_base-Lz))/Lz;
// 		
// 		plan_fftw();
// 		make_C();
// 		this->set_IC();
// 		this->evolve();
// 		ED_vals_Lz(count) = ED;		
// 		count++;
// 	}
// 	
// 	Lz = Lz_base;
// 	dz = dz_base;
// 	count = 0;
	
// 	for (int i=dm; i<=dp; i++)
// 	{
// 		Lx = Lx_base + inc*Lx_base*i;
// 		Ly = Ly_base + inc*Ly_base*i;
// 		Lz = Lz_base + inc*Lz_base*i;
// 		
// 		dx = Lx/(Nx);
// 		dy = Ly/(Ny);
// 		dz = Lz/(Nz);
// 		
// 		double V_base = Lx_base*Ly_base*Lz_base; 
// 		double V = Lx*Ly*Lz;
// 		
// 	 	//n_bar = (n_base*V_base + (V_base-V))/V;
// 	 	
// 		plan_fftw();
// 		make_C();
// 		this->set_IC();
// 		cout << "n_bar: " << n_bar << endl;
// 		//psi1 = psi1+n1;
// 		//test_mass = sum(real(psi1))*dx*dy*dz;
// 		test_mass = sum(real(psi1)+1.0)*dx*dy*dz;
//  		cout << "test mass: " << test_mass << endl;
//  		psi_mean = mean(real(psi1));
//  		cout << "psi mean: " << psi_mean << endl;
// 		this->evolve();
// 		ED_vals_V(count) = ED;
// 		count++;
// 		write_out(real(psi1),"psi" + to_string(count));
// 			
// 	}
	int i_stop = dp-dm;
	cout << "First val: " << ED_vals_V(0) << " Last val: " << ED_vals_V(i_stop) << " Min val: " << min(ED_vals_V) <<endl;
	ofstream strainFile;
 	strainFile.open("strainData.txt");
 	strainFile.precision(16);
 	strainFile.setf(std::ios::fixed, std:: ios::floatfield); // floatfield set to fixed
//   	
//   	
//   	//strainFile << "Line info: (1)Lx_strain, (2)Ly_strain, (3)Lz_strain, (4)Lx_start, Lx_inc, Lx_stop, (5)Ly_start, Ly_inc, Ly_stop, (6)Lz_start, Lz_inc, Lz_stop " << endl << endl;
  	for (int i = 0; i<i_stop; i++) { strainFile << ED_vals_Lx(i) << " ";}
  	strainFile << endl;
//   	
//   	//strainFile << "Ly_strain:" << endl;
//   	for (int i = 0; i<i_stop; i++) { strainFile << ED_vals_Ly(i) << " ";}
//   	strainFile << endl;
//   	
//   	//strainFile << "Lz_strain:" << endl;
//   	for (int i = 0; i<i_stop; i++) { strainFile << ED_vals_Lz(i) << " ";}
//   	strainFile << endl;
//   	
//   	//strainFile << "V_strain:" << endl;
//   	for (int i = 0; i<=i_stop; i++) { strainFile << ED_vals_V(i) << " ";}
//   	strainFile << endl << endl;
//   
//   	strainFile << Lx_base + dm*inc*Lx_base << " " << Lx_base*inc << " " << Lx_base + dp*inc*Lx_base << endl;
// 	strainFile << Ly_base + dm*inc*Ly_base << " " << Ly_base*inc << " " << Ly_base + dp*inc*Ly_base << endl;
// 	strainFile << Lz_base + dm*inc*Lz_base << " " << Lz_base*inc << " " << Lz_base + dp*inc*Lz_base << endl;
//   
 	strainFile.close();
// 	
// 	cout << ED_vals_Lx << ED_vals_Ly << ED_vals_Lz << ED_vals_V << endl;
// 	cout << "Lx start: " << Lx_base + dm*inc*Lx_base << " x_inc: " << Lx_base*inc << " Lx stop: " <<  Lx_base + dp*inc*Lx_base << endl;
// 	cout << "Ly start: " << Ly_base + dm*inc*Ly_base << " y_inc: " << Ly_base*inc << " Ly stop: " <<  Ly_base + dp*inc*Ly_base << endl;
// 	cout << "Lz start: " << Lz_base + dm*inc*Lz_base << " z_inc: " << Lz_base*inc << " Lz stop: " <<  Lz_base + dp*inc*Lz_base << endl<< endl;
}