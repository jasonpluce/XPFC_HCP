#include "header.h"
#include <blitz/array.h>
#include <math.h>
#include <time.h>

using namespace blitz;
using namespace std;

void GB_system::set_IC() 
{
	cout << "Creating initial conditions" << endl;
	srand(time(NULL));
	
	Array<double,1> x_grid(Nx);
	Array<double,1> y_grid(Ny);
	Array<double,1> z_grid(Nz);
	
	double qt = q;
	double x_val, x_val2;
	double y_val, y_val2;
	double z_val, z_val2;
	double psi_mean = 0.0;
	double temp1,temp2;
	
 	firstIndex i1;
 	x_grid = (i1)*dx;
 	y_grid = (i1)*dy; 
 	z_grid = (i1)*dz;
	
	//cout << x_grid << y_grid << z_grid;
	
	int offset = IC_spacing;
	
	if (IC_type == "HEXAGONAL"){
		for (int i=0; i< Nz; i++){	
			for (int j=(0+offset);j<(Ny/2-offset);j++){
				for (int k = (0+offset);k<(Nx-offset);k++){	
				
					//cout << "i: " << i << " j: " << j << " k: " << k << endl;
					x_val = x_grid(k)*cos(theta) - y_grid(j)*sin(theta);
					y_val = x_grid(k)*sin(theta) + y_grid(j)*cos(theta);
					
					
					psi1(i,j,k) = (0.5*(cos(2*qt*y_val/sqrt(3))) + cos(qt*x_val)*cos(qt*y_val/sqrt(3))); // OMA of 3D hexagonal density profile
				}
 			}
 		}

 		theta = -theta;
 		
		for (int i=0; i< Nz; i++){	
			for (int j=(Ny/2+offset); j<(Ny-offset); j++){
				for (int k = (0+offset);k<(Nx-offset);k++){	
				
					x_val = x_grid(k)*cos(theta) - y_grid(j)*sin(theta);
					y_val = x_grid(k)*sin(theta) + y_grid(j)*cos(theta);
					
					psi1(i,j,k) = (0.5*(cos(2*qt*y_val/sqrt(3))) + cos(qt*x_val)*cos(qt*y_val/sqrt(3))); // OMA of 3D hexagonal density profile
				}
 			}
		}
 	}
	

	else if (IC_type == "HCP_ZROT"){
	
		for (int i=(0+offset);i<(Nz-offset);i++){
			for (int j=(0+offset);j<(Ny/2-offset);j++){
				for (int k=(0+offset); k<(Nx-offset); k++){
					
					x_val = x_grid(k)*cos(theta) - y_grid(j)*sin(theta);
					y_val = x_grid(k)*sin(theta) + y_grid(j)*cos(theta);
					z_val = z_grid(i);
					
					x_val2 = x_val;
					y_val2 = y_val + 1/sqrt(3.0);
					z_val2 = z_val + 0.5*c;
			
					temp1 = (0.5*(cos(2*qt*y_val/sqrt(3.0))) + cos(qt*x_val)*cos(qt*y_val/sqrt(3.0)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val)); // OMA of 3D hexagonal density profile
					temp2 = (0.5*(cos(2*qt*y_val2/sqrt(3.0))) + cos(qt*x_val2)*cos(qt*y_val2/sqrt(3.0)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val2));
					psi1(i,j,k) = max(temp1,temp2);
					//psi1(i,j,k) = (double) rand()/(RAND_MAX);
				}
 			}
 		}
 
 		theta = -theta;
 	
 		for (int i=(0+offset); i<(Nz-offset); i++){
			for (int j=(Ny/2+offset); j<(Ny-offset); j++){
				for (int k=(0+offset); k<(Nx-offset); k++){
		
					x_val = x_grid(k)*cos(theta) - y_grid(j)*sin(theta);
					y_val = x_grid(k)*sin(theta) + y_grid(j)*cos(theta);
					z_val = z_grid(i);
					
					x_val2 = x_val;
					y_val2 = y_val + 1/sqrt(3.0);
					z_val2 = z_val + 0.5*c;
			
					
					temp1 = (0.5*(cos(2*qt*y_val/sqrt(3.0))) + cos(qt*x_val)*cos(qt*y_val/sqrt(3.0)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val)); // OMA of 3D hexagonal density profile
					temp2 = (0.5*(cos(2*qt*y_val2/sqrt(3.0))) + cos(qt*x_val2)*cos(qt*y_val2/sqrt(3.0)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val2));
					
					psi1(i,j,k) = max(temp1,temp2);
					//psi1(i,j,k) = (double) rand()/(RAND_MAX);
				}
 			}
		}
	}
	else if (IC_type == "HCP_YROT"){
		int x_shift, y_shift, z_shift;
		// D_term = 0.25*sqrt(3.0);
		for (int i=(0+offset);i<(Nz-offset);i++){	
			for (int j=(0+offset);j<(Ny-offset);j++){
				for (int k=(0+offset); k<(Nx/2-offset); k++){
								
					x_val = x_grid(k)*cos(theta) - z_grid(i)*sin(theta);					
					y_val = y_grid(j);
					z_val = x_grid(k)*sin(theta) + z_grid(i)*cos(theta);
					
					x_val2 = x_val;
					y_val2 = y_val + 1.0/sqrt(3.0);
					z_val2 = z_val + 0.5*c;
								
					temp1 = (0.5*(cos(2*qt*y_val/sqrt(3.0))) + cos(qt*x_val)*cos(qt*y_val/sqrt(3.0)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val)); // OMA of 3D hexagonal density profile								
					temp2 = (0.5*(cos(2*qt*y_val2/sqrt(3.0))) + cos(qt*x_val2)*cos(qt*y_val2/sqrt(3.0)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val2));
					
					//Alternative IC (does not work)
					//temp1 = cos(qt*x_val) + 2*cos(0.5*qt*y_val*sqrt(3.0))*cos(0.5*qt*x_val) + cos(qt*sqrt(3.0/8.0)*z_val); // OMA of 3D hexagonal density profile
					//temp2 = cos(qt*x_val2) + 2*cos(0.5*qt*y_val2*sqrt(3.0))*cos(0.5*qt*x_val2) + cos(qt*sqrt(3.0/8.0)*z_val2);

					psi1(i,j,k) = max(temp1,temp2);
					//psi1(i,j,k) = temp1 + temp2;
					//psi1(i,j,k) = (double) rand()/(RAND_MAX);
				}
 			}
 		}

 		theta = -theta;
 		
 		x_shift = 1.0;
 		y_shift = 0.0;
 		z_shift = 0.0;
 		
 		for (int i=(0+offset); i<(Nz-offset); i++){
			for (int j=(0+offset); j<(Ny-offset); j++){
				for (int k=(Nx/2+offset); k<(Nx-offset); k++){
		
					x_val = x_grid(k + x_shift)*cos(theta) - z_grid(i + z_shift)*sin(theta);
					y_val = y_grid(j + y_shift);
					z_val = x_grid(k + x_shift)*sin(theta) + z_grid(i + z_shift)*cos(theta);

					x_val2 = x_val;				
					y_val2 = y_grid(j) + 1.0/sqrt(3.0);
					z_val2 = z_val + 0.5*c;
			
					temp1 = (0.5*(cos(2*qt*y_val/sqrt(3.0))) + cos(qt*x_val)*cos(qt*y_val/sqrt(3.0)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val)); // OMA of 3D hexagonal density profile
					temp2 = (0.5*(cos(2*qt*y_val2/sqrt(3.0))) + cos(qt*x_val2)*cos(qt*y_val2/sqrt(3.0)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val2));
					
					//Alternative IC (does not work)
					//temp1 = cos(qt*x_val) + 2*cos(0.5*qt*y_val*sqrt(3.0))*cos(0.5*qt*x_val) + cos(qt*sqrt(3.0/8.0)*z_val); // OMA of 3D hexagonal density profile
					//temp2 = cos(qt*x_val2) + 2*cos(0.5*qt*y_val2*sqrt(3.0))*cos(0.5*qt*x_val2) + cos(qt*sqrt(3.0/8.0)*z_val2);

					psi1(i,j,k) = max(temp1,temp2);
					//psi1(i,j,k) = (double) rand()/(RAND_MAX);
					
				}
 			}
		}
	}
	else if (IC_type == "HCP_UC"){
	
	 	//offset = 0.1*Nx;
		//cout << "Offset 0.1*Nx: " << offset;
	 	
		offset = 5.0;
	 	//target_angle = 0.0;
		//theta = 0.0;
	
		for (int i=(0+offset);i<(Nz-offset);i++){
			for (int j=(0+offset);j<(Ny-offset);j++){
				for (int k=(0+offset); k<(Nx-offset); k++){
					
					x_val = x_grid(k)*cos(theta) - y_grid(j)*sin(theta);
					y_val = x_grid(k)*sin(theta) + y_grid(j)*cos(theta);
					z_val = z_grid(i);
					
					x_val2 = x_val;
					y_val2 = y_val + 1/sqrt(3.0);
					z_val2 = z_val + 0.5*c;
					
					temp1 = (0.5*(cos(2*qt*y_val/sqrt(3))) + cos(qt*x_val)*cos(qt*y_val/sqrt(3)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val)); // OMA of 3D hexagonal density profile
					temp2 = (0.5*(cos(2*qt*y_val2/sqrt(3))) + cos(qt*x_val2)*cos(qt*y_val2/sqrt(3)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val2));
					
					psi1(i,j,k) = max(temp1,temp2);
					//psi1(i,j,k) = temp1;
					//psi1(i,j,k) = (double) rand()/(RAND_MAX);
				}
 			}
 		}
	}
 	else if (IC_type == "BCC_UC"){
 		
 		offset = 3.0;
 		target_angle = 0.0;
		//theta = DEG2RAD*45.0;
		qt = 2.0*PI/a;
 		
 		for (int i=(0+offset);i<(Nz-offset);i++){
			for (int j=(0+offset);j<(Ny-offset);j++){
				for (int k=(0+offset); k<(Nx-offset); k++){
		
					x_val = x_grid(k)*cos(theta) - y_grid(j)*sin(theta);
					y_val = x_grid(k)*sin(theta) + y_grid(j)*cos(theta);
					z_val = z_grid(i);
									
					psi1(i,j,k) = (cos(qt*y_val)*cos(qt*x_val)+cos(qt*y_val)*cos(qt*z_val)+cos(qt*x_val)*cos(qt*z_val));
					psi1(i,j,k) = (double) rand()/(RAND_MAX);
				}
 			}
 		}
 	}
 	
 	else if (IC_type == "FCC_UC"){
 		
 		offset = 3.0;
 		target_angle = 0.0;
		theta = 0.0;
		qt = 2.0*PI/a;
 		
 		for (int i=(0+offset);i<(Nz-offset);i++){
			for (int j=(0+offset);j<(Ny-offset);j++){
				for (int k=(0+offset); k<(Nx-offset); k++){
		
					x_val = x_grid(k)*cos(theta) - y_grid(j)*sin(theta);
					y_val = x_grid(k)*sin(theta) + y_grid(j)*cos(theta);
					z_val = z_grid(i);
					
					psi1(i,j,k) = cos(qt*x_val)*cos(qt*y_val)*cos(qt*z_val);// + cos(0.5*qt*x_val)*cos(0.5*qt*y_val)*cos(0.5*qt*z_val);
					//psi1(i,j,k) = (cos(qt*y_val)*cos(qt*x_val)*cos(qt*z_val)+ cos(2.0*qt*y_val) + cos(2.0*qt*x_val) + cos(2.0*qt*z_val));
					//psi1(i,j,k) = (double) rand()/(RAND_MAX);
				}
 			}
 		}
 	}
 	
 	 else if (IC_type == "SC_UC"){
 		
 		offset = 0;
 		target_angle = 0.0;
		theta = 0.0;
 		
 		for (int i=(0+offset);i<(Nz-offset);i++){
			for (int j=(0+offset);j<(Ny-offset);j++){
				for (int k=(0+offset); k<(Nx-offset); k++){
		
					x_val = x_grid(k)*cos(theta) - y_grid(j)*sin(theta);
					y_val = x_grid(k)*sin(theta) + y_grid(j)*cos(theta);
					z_val = z_grid(i);
	
					psi1(i,j,k) = cos(qt*x_val)*cos(qt*y_val)*cos(qt*z_val);
					//psi1(i,j,k) = (double) rand()/(RAND_MAX);
				}
 			}
 		}
 	}
 	
 	else if (IC_type == "Columnar"){
 		
 		offset = 0.0;
 		double side_angle = 0.0;
 		double middle_angle = 4.0;
 		
 		target_angle = 0.0;
		double theta1 = DEG2RAD*side_angle;
		double theta2 = -DEG2RAD*side_angle;
		double theta3 = DEG2RAD*middle_angle;
		//theta2 = 0.0;
 		double R = 0.4*min(Nx,Ny);
 		double hNx = 0.5*Nx;
 		double hNy = 0.5*Ny;
 		
 		for (int i=(0+offset);i<(Nz-offset);i++){
			for (int j=(0+offset);j<(Ny-offset);j++){
				for (int k=(0+offset); k<(Nx-offset); k++){
					
					//if ((j > hNy-R) && (j < hNy+R) && (k > hNx-R) && (k < hNx+R))

					
					if (sqrt((hNy-j)*(hNy-j) + (hNx-k)*(hNx-k)) > R + offset)
 					{
						if (k > hNx + offset){
							x_val = x_grid(k)*cos(theta1) - y_grid(j)*sin(theta1);
							y_val = x_grid(k)*sin(theta1) + y_grid(j)*cos(theta1);
							z_val = z_grid(i);
					
							x_val2 = x_val;
							y_val2 = y_val + 1/sqrt(3.0);
							z_val2 = z_val + 0.5*c;
						
							temp1 = (0.5*(cos(2*qt*y_val/sqrt(3))) + cos(qt*x_val)*cos(qt*y_val/sqrt(3)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val)); // OMA of 3D hexagonal density profile
							temp2 = (0.5*(cos(2*qt*y_val2/sqrt(3))) + cos(qt*x_val2)*cos(qt*y_val2/sqrt(3)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val2));
							psi1(i,j,k) = max(temp1,temp2);
							
							//psi1(i,j,k) = (cos(qt*y_val)*cos(qt*x_val)+cos(qt*y_val)*cos(qt*z_val)+cos(qt*x_val)*cos(qt*z_val));//BCC
						}
						
						
						else if (k <= hNx - offset)
						{
							x_val = x_grid(k)*cos(theta2) - y_grid(j)*sin(theta2);
							y_val = x_grid(k)*sin(theta2) + y_grid(j)*cos(theta2);
							z_val = z_grid(i);
					
							x_val2 = x_val;
							y_val2 = y_val + 1/sqrt(3.0);
							z_val2 = z_val + 0.5*c;
					
							temp1 = (0.5*(cos(2*qt*y_val/sqrt(3))) + cos(qt*x_val)*cos(qt*y_val/sqrt(3)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val)); // OMA of 3D hexagonal density profile
							temp2 = (0.5*(cos(2*qt*y_val2/sqrt(3))) + cos(qt*x_val2)*cos(qt*y_val2/sqrt(3)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val2));
							psi1(i,j,k) = max(temp1,temp2);
							
							//psi1(i,j,k) = (cos(qt*y_val)*cos(qt*x_val)+cos(qt*y_val)*cos(qt*z_val)+cos(qt*x_val)*cos(qt*z_val));// BCC
						}
					}
						
					else if (sqrt((hNy-j)*(hNy-j) + (hNx-k)*(hNx-k)) < R - offset)
					
					{
						x_val = x_grid(k)*cos(theta3) - y_grid(j)*sin(theta3);
						y_val = x_grid(k)*sin(theta3) + y_grid(j)*cos(theta3);
						z_val = z_grid(i);
					
						x_val2 = x_val;
						y_val2 = y_val + 1/sqrt(3.0);
						z_val2 = z_val + 0.5*c;
					
						temp1 = (0.5*(cos(2*qt*y_val/sqrt(3))) + cos(qt*x_val)*cos(qt*y_val/sqrt(3)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val)); // OMA of 3D hexagonal density profile
						temp2 = (0.5*(cos(2*qt*y_val2/sqrt(3))) + cos(qt*x_val2)*cos(qt*y_val2/sqrt(3)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val2));
						psi1(i,j,k) = max(temp1,temp2);
						
						//psi1(i,j,k) = (cos(qt*y_val)*cos(qt*x_val)+cos(qt*y_val)*cos(qt*z_val)+cos(qt*x_val)*cos(qt*z_val));//BCC
					}	
					//psi1(i,j,k) = (double) rand()/(RAND_MAX);
				}
 			}
 		}
 	}
 	
 	else if (IC_type == "Spherical"){
 		
 		offset = 0.0;
 		double side_angle = 0.0;
 		double middle_angle = 20.0;
 		
 		target_angle = 0.0;
		double theta1 = DEG2RAD*side_angle;
		double theta2 = -DEG2RAD*side_angle;
		double theta3 = DEG2RAD*middle_angle;
		//theta2 = 0.0;
 		double R = 0.4*min(min(Nx,Ny),Nz);
 		double hNx = 0.5*Nx;
 		double hNy = 0.5*Ny;
 		double hNz = 0.5*Nz;
 		
 		double ux = 1.0/sqrt(3.0);
 		double uy = 1.0/sqrt(3.0);
 		double uz = 1.0/sqrt(3.0);
 		
 		for (int i=(0+offset);i<(Nz-offset);i++){
			for (int j=(0+offset);j<(Ny-offset);j++){
				for (int k=(0+offset); k<(Nx-offset); k++){
										
					if (sqrt((hNz-i)*(hNz-i) + (hNy-j)*(hNy-j) + (hNx-k)*(hNx-k)) > R + offset)
 					{
						x_val = x_grid(k)*cos(theta1) - y_grid(j)*sin(theta1);
						y_val = x_grid(k)*sin(theta1) + y_grid(j)*cos(theta1);
						z_val = z_grid(i);
												
						x_val2 = x_val;
						y_val2 = y_val + 1/sqrt(3.0);
						z_val2 = z_val + 0.5*c;
						
						temp1 = (0.5*(cos(2*qt*y_val/sqrt(3))) + cos(qt*x_val)*cos(qt*y_val/sqrt(3)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val)); // OMA of 3D hexagonal density profile
						temp2 = (0.5*(cos(2*qt*y_val2/sqrt(3))) + cos(qt*x_val2)*cos(qt*y_val2/sqrt(3)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val2));
						psi1(i,j,k) = max(temp1,temp2);
							
						//psi1(i,j,k) = (cos(qt*y_val)*cos(qt*x_val)+cos(qt*y_val)*cos(qt*z_val)+cos(qt*x_val)*cos(qt*z_val));//BCC
					}	
					else
					{
// 						x_val = x_grid(k)*cos(theta3) - y_grid(j)*sin(theta3);
// 						y_val = x_grid(k)*sin(theta3) + y_grid(j)*cos(theta3);
// 						z_val = z_grid(i);
					
						x_val = x_grid(k)*(cos(theta3) + ux*ux*(1-cos(theta3))) + y_grid(j)*(ux*uy*(1-cos(theta3)) - uz*sin(theta3)) + z_grid(i)*(ux*uz*(1.0-cos(theta3)) + uy*sin(theta3));
						y_val = x_grid(k)*(uy*ux*(1.0-cos(theta3)) + uz*sin(theta3)) + y_grid(j)*(cos(theta3) + uy*uy*(1.0-cos(theta3))) + z_grid(i)*(uy*uz*(1.0-cos(theta3)) - ux*sin(theta3));
						z_val = x_grid(k)*(uz*ux*(1.0-cos(theta3)) - uy*sin(theta3)) - y_grid(j)*(uy*uz*(1.0-cos(theta3)) + ux*sin(theta3)) + z_grid(i)*(cos(theta3) + uz*uz*(1-cos(theta3)));
					
						x_val2 = x_val;
						y_val2 = y_val + 1/sqrt(3.0);
						z_val2 = z_val + 0.5*c;
					
						temp1 = (0.5*(cos(2*qt*y_val/sqrt(3))) + cos(qt*x_val)*cos(qt*y_val/sqrt(3)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val)); // OMA of 3D hexagonal density profile
						temp2 = (0.5*(cos(2*qt*y_val2/sqrt(3))) + cos(qt*x_val2)*cos(qt*y_val2/sqrt(3)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val2));
						psi1(i,j,k) = max(temp1,temp2);
						
						//psi1(i,j,k) = (cos(qt*y_val)*cos(qt*x_val)+cos(qt*y_val)*cos(qt*z_val)+cos(qt*x_val)*cos(qt*z_val));//BCC
					}
					
					//psi1(i,j,k) = (double) rand()/(RAND_MAX);
				}
			}
		}
 	}
 						//if ((j > hNy-R) && (j < hNy+R) && (k > hNx-R) && (k < hNx+R))
// 					if (sqrt((hNz-i)*(hNz-i) + (hNy-j)*(hNy-j) + (hNx-k)*(hNx-k)) > R + offset)
//  					{
// 						if (k > hNx + offset){
// 							x_val = x_grid(k)*cos(theta1) - y_grid(j)*sin(theta1);
// 							y_val = x_grid(k)*sin(theta1) + y_grid(j)*cos(theta1);
// 							z_val = z_grid(i);
// 												
// 							x_val2 = x_val;
// 							y_val2 = y_val + 1/sqrt(3.0);
// 							z_val2 = z_val + 0.5*c;
// 						
// 							temp1 = (0.5*(cos(2*qt*y_val/sqrt(3))) + cos(qt*x_val)*cos(qt*y_val/sqrt(3)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val)); // OMA of 3D hexagonal density profile
// 							temp2 = (0.5*(cos(2*qt*y_val2/sqrt(3))) + cos(qt*x_val2)*cos(qt*y_val2/sqrt(3)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val2));
// 							psi1(i,j,k) = max(temp1,temp2);
// 							
// 							//psi1(i,j,k) = (cos(qt*y_val)*cos(qt*x_val)+cos(qt*y_val)*cos(qt*z_val)+cos(qt*x_val)*cos(qt*z_val));//BCC
// 						}
// 						
// 						
// 						else if (k <= hNx - offset)
// 						{
// 							x_val = x_grid(k)*cos(theta2) - y_grid(j)*sin(theta2);
// 							y_val = x_grid(k)*sin(theta2) + y_grid(j)*cos(theta2);
// 							z_val = z_grid(i);
// 					
// 							x_val2 = x_val;
// 							y_val2 = y_val + 1/sqrt(3.0);
// 							z_val2 = z_val + 0.5*c;
// 					
// 							temp1 = (0.5*(cos(2*qt*y_val/sqrt(3))) + cos(qt*x_val)*cos(qt*y_val/sqrt(3)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val)); // OMA of 3D hexagonal density profile
// 							temp2 = (0.5*(cos(2*qt*y_val2/sqrt(3))) + cos(qt*x_val2)*cos(qt*y_val2/sqrt(3)) + 2.0*cos(qt*sqrt(3.0/8.0)*z_val2));
// 							psi1(i,j,k) = max(temp1,temp2);
// 							
// 							//psi1(i,j,k) = (cos(qt*y_val)*cos(qt*x_val)+cos(qt*y_val)*cos(qt*z_val)+cos(qt*x_val)*cos(qt*z_val));// BCC
// 						}
// 					}
 	
 	 else if (IC_type == "Liquid"){
 		
 		for (int i=(0+offset);i<(Nz-offset);i++){
			for (int j=(0+offset);j<(Ny-offset);j++){
				for (int k=(0+offset); k<(Nx-offset); k++){
		
					psi1(i,j,k) = n_bar;
					//psi1(i,j,k) = (double) rand()/(RAND_MAX);
				}
 			}
 		}
 	}

 	double pmin = min(real(psi1));
 	double pmax = max(real(psi1));
 	
 	if (IC_type != "Liquid"){psi1 = 2.0*(psi1-pmin)/(pmax-pmin) - 0.5;}
 	
 	psi_mean = mean(real(psi1)); // determine average density of the system
 	//cout << "psi_mean: " << psi_mean << endl;
 	psi1 = psi1 - psi_mean + n_bar; // set average value of psi = n_bar
 	cout << "psi_mean2: " << mean(real(psi1)) << endl; 
 	// mass = sum(real(psi1))*dx*dy*dz;// determine mass of system
    // cout << "mass: " << mass << endl;
 	write_out(real(psi1),"IC");
 	cout << endl;
 	
}
