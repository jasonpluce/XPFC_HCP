#include "header.h"
#include <blitz/array.h>
#include <math.h>
#include <fstream>

void GB_system::make_C()
{

	Array<double,3> k_vals(Nz,Ny,Nx);
	
	Array<double,1> k1x(Nx);
	Array<double,1> k1y(Ny);
	Array<double,1> k1z(Nz);
	
	Array<int,1> kx(Nx);
	Array<int,1> ky(Ny);
	Array<int,1> kz(Nz);
	
	double paramx = 2.0*PI/Lx;
	double paramy = 2.0*PI/Ly;
	double paramz = 2.0*PI/Lz;
	
	double k0 = 0.0;
	
	double temp0, temp1, temp2, temp3, temp4;
	
	for (int i=0; i<=Nx/2; i++){
		k1x(i) = i*paramx;
		kx(i) = i;
	}

	for (int i=Nx/2; i>0; i--){
		k1x(Nx-i) = (-i)*paramx;
		kx(Nx-i) = -i;
	}

	for (int j=0; j<=Ny/2; j++){
		k1y(j) = j*paramy;
		ky(j) = j;
	}

	for (int j=Ny/2; j>0; j--){
		k1y(Ny-j) = (-j)*paramy;
		ky(Ny-j) = -j;
	}
	
	for (int k=0; k<=Nz/2; k++){
		k1z(k) = k*paramz;
		kz(k) = k;
	}

	for (int k=Nz/2; k>0; k--){
		k1z(Nz-k) = (-k)*paramz;
		kz(Nz-k) = -k;
	}
	
	for (int i=0; i<Nz; i++)
	{
    	for (int j=0; j<Ny; j++)
    	{
    		for(int k=0; k<Nx; k++)
    		{
    			temp1 = pow2(k1x(k)) + pow2(k1y(j)) + pow2(k1z(i));
    			ks(i,j,k) = temp1;
    			k_vals(i,j,k) = sqrt(temp1);
    		}
 		}
 	}
	
	//temp1 = 0;
 	//write_out(ks,"ks.txt");
 	//cout << k_vals.size() << endl;
 	//cout << k_vals << endl;
 	//write_out(k_vals,"k_vals");
 	
	if (corr_type == "HEXAGONAL"){
		cout << "Generating HEXAGONAL correlation function" << endl;
		firstIndex i;
    	secondIndex j;
    	thirdIndex l;
		
		
 		k1 = 4.0*PI/(a*sqrt(3.0));
 		k2 = 4.0*PI;
 		k3 = 8.0*PI/(a*sqrt(3.0));
 		
 		C = exp(-pow2(k_vals(i,j,l)-k1)/(2.0*(alpha1*alpha1)))*exp(-(sigma*sigma*k1*k1)/(2.0*rho1*beta1));
 		
// 		for (int i=0; i<Nz; i++)
// 		{
//     		for (int j=0; j<Ny; j++)
//     		{
//     			for(int k=0; k<Nx; k++)
//     			{
//      			temp0 = exp(-pow2(k_vals(i,j,k)-k0)/(2.0*(alpha1*alpha1)));
//      			temp1 = exp(-pow2(k_vals(i,j,k)-k1)/(2.0*(alpha1*alpha1)))*exp(-(sigma*sigma*k1*k1)/(2.0*rho1*beta1));
//      			temp2 = exp(-pow2(k_vals(i,j,k)-k2)/(2.0*(alpha1*alpha1)))*exp(-(sigma*sigma*k2*k2)/(2.0*rho1*beta1));
//      			temp3 = exp(-pow2(k_vals(i,j,k)-k3)/(2.0*(alpha1*alpha1)))*exp(-(sigma*sigma*k3*k3)/(2.0*rho1*beta1));
// 					temp2 = 0.0;
// 					temp3 = 0.0;
// 					//temp0 = 40.0*temp0;
// 					temp0 = 0.0;
//     				C(i,j,k) = max(max(temp1,temp2), temp3) - temp0;
//     			}
//  		}
//  	}
	}
	
	else if (corr_type == "HCP"){
	
    	cout << "Generating HCP correlation function" << endl;
//		ofstream HCP_DCF;
// 		ofstream HCP_DCF_P1;
// 		ofstream HCP_DCF_P2;
// 		ofstream HCP_DCF_P3;
// 		ofstream HCP_DCF_P0;
// 		
// 		HCP_DCF.open("HCP_DCF.txt", ofstream::app);
// 		HCP_DCF.precision(16);
// 		HCP_DCF.setf(std::ios::fixed, std:: ios::floatfield); // floatfield set to fixed
// 		
// 		HCP_DCF_P1.open("HCP_DCF_P1.txt", ofstream::app);
// 		HCP_DCF_P1.precision(16);
// 		HCP_DCF_P1.setf(std::ios::fixed, std:: ios::floatfield); // floatfield set to fixed
// 		
// 		HCP_DCF_P2.open("HCP_DCF_P2.txt", ofstream::app);
// 		HCP_DCF_P2.precision(16);
// 		HCP_DCF_P2.setf(std::ios::fixed, std:: ios::floatfield); // floatfield set to fixed
// 		
// 		HCP_DCF_P3.open("HCP_DCF_P3.txt", ofstream::app);
// 		HCP_DCF_P3.precision(16);
// 		HCP_DCF_P3.setf(std::ios::fixed, std:: ios::floatfield); // floatfield set to fixed
// 		
// 		HCP_DCF_P0.open("HCP_DCF_P0.txt", ofstream::app);
// 		HCP_DCF_P0.precision(16);
// 		HCP_DCF_P0.setf(std::ios::fixed, std:: ios::floatfield); // floatfield set to fixed
		
		double k1_scale = 1.0;
		double k2_scale = 1.0;
		double k3_scale = 1.0;
		
 		k1 = 4.0*PI/(sqrt(3.0));
 		k2 = 4.0*PI*sqrt(3.0/8.0);
 		k3 = 2.0*PI*sqrt(41.0)/(sqrt(24.0));
		//double k4 = 4.0*PI/a;
    	for (int i=0; i<Nz; i++)
		{
    		for (int j=0; j<Ny; j++)
    		{
    			for(int k=0; k<Nx; k++)
    			{
    				//temp0 = exp(-pow2(k_vals(i,j,k)-k0)/(2.0*(alpha1*alpha1)));
     				temp1 = k1_scale*exp(-pow2(k_vals(i,j,k)-k1)/(2.0*(alpha1*alpha1)))*exp(-(sigma*sigma*k1*k1)/(2.0*rho1*beta1));
     				temp2 = k2_scale*exp(-pow2(k_vals(i,j,k)-k2)/(2.0*(alpha1*alpha1)))*exp(-(sigma*sigma*k2*k2)/(2.0*rho2*beta2));
     				temp3 = k3_scale*exp(-pow2(k_vals(i,j,k)-k3)/(2.0*(alpha1*alpha1)))*exp(-(sigma*sigma*k3*k3)/(2.0*rho3*beta3));
					//temp4 = exp(-pow2(k_vals(i,j,k)-k4)/(2.0*(alpha1*alpha1)))*exp(-(sigma*sigma*k4*k4)/(2.0*rho3*beta3));
					//temp0 = -40.0*temp0;
					//temp0 = 0.0;
					//C(i,j,k) = temp1;
					//C(i,j,k) = temp2;
					//C(i,j,k) = temp3;
					//C(i,j,k) = max(temp2,temp3);
					//C(i,j,k) = max(temp1,temp3);
					//C(i,j,k) = max(temp1,temp2);// + temp0;
					C(i,j,k) = max(max(temp1,temp2), temp3);// + temp0;
    				//C(i,j,k) = temp3;
    				//C(i,j,k) = max(temp2,temp4);// + temp0;
					//C(i,j,k) = temp4;
   					//HCP_DCF << kz(i) << "\t" << ky(j) << "\t" << kx(k) << "\t" << real(C(i,j,k)) << endl;
//   					HCP_DCF_P1 << kz(i) << "\t" << ky(j) << "\t" << kx(k) << "\t" << temp1 << endl;
//   					HCP_DCF_P2 << kz(i) << "\t" << ky(j) << "\t" << kx(k) << "\t" << temp2 << endl;
//   					HCP_DCF_P3 << kz(i) << "\t" << ky(j) << "\t" << kx(k) << "\t" << temp3 << endl;
//   					HCP_DCF_P0 << kz(i) << "\t" << ky(j) << "\t" << kx(k) << "\t" << temp0 << endl;
    			}
 			}
 		}
  		//HCP_DCF.close();
//  		HCP_DCF_P1.close();
//  		HCP_DCF_P2.close();
//  		HCP_DCF_P3.close();
//  		HCP_DCF_P0.close();
 	}
 	
 	else if (corr_type == "BCC"){
 	
 	    cout << "Generating BCC correlation function" << endl;
 		k1 = 2.0*PI*sqrt(2.0)/a;
 		k2 = 4.0*PI/a;
 		alpha1 = 1.0;
 		alpha2 = 1.0;
 		beta1 = 12.0;
 		beta2 = 6.0;
 		sigma = 0.12;
 		n_bar = 0.15;
 		
 		for (int i=0; i<Nz; i++)
		{
    		for (int j=0; j<Ny; j++)
    		{
    			for(int k=0; k<Nx; k++)
    			{
    				temp1 = exp(-pow2(k_vals(i,j,k)-k1)/(2.0*(alpha1*alpha1)))*exp(-(sigma*sigma*k1*k1)/(2.0*rho1*beta1));
     				temp2 = exp(-pow2(k_vals(i,j,k)-k2)/(2.0*(alpha2*alpha2)))*exp(-(sigma*sigma*k2*k2)/(2.0*rho2*beta2));
     				C(i,j,k) = max(temp1,temp2);
    			}
 			}
 		}	
 	}
 	else
 	{
 		cout << "Undefined correlation function";
 		exit (EXIT_FAILURE);
 	}
 	//cout << "Writing correlation function to \"C.vtk\"" << endl;
  	write_out(real(C),"C"); 
} 
