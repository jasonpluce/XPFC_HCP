//#include "functions.h"
#include "header.h"
#include <blitz/array.h>
#include <math.h>
#include <iostream>

using namespace blitz;
using namespace std;

void GB_system::get_size() 
{
	cout << "Determining system size" << endl;
	
	Array<double,1> a_vect(2), b_vect(2), temp1(2), temp2(2);

	double target_theta = target_angle*DEG2RAD;// convert target angle to radians
	cout << "target angle: " << target_angle << endl;
    double theta2;
    double da = 1e15;
    
	double i = 0;
	double j = 0;
	double m = 0;
	double l = 0;
	
	double l1 = 1.0;
	double l3 = 1.0;
	double m1 = 2.0;
	double n1 = 1.0;
	double n3 = 1.0;
	
	double c_num = 1.0;
	
	//double temp_Lx;
	//double temp_Ly;
	//double temp_Lz;
	
	//int a_num = 3.0; // number of unit cells in the z-direction
	
	
	
	if (IC_type == "HEXAGONAL"){
		cout << "Caluculating hexagonal 2D system size (Lx, Ly), rotated about z-axis" << endl; 
		
		
		l = 5;
		while ((da > da_thresh) && (l<=50))
		{
			l = l+1;
			m = round(2.0*l*tan(target_theta)/(tan(target_theta)-sqrt(3.0)));
			theta = atan((sqrt(3.0)*m)/(m-2*l));
			da = abs(RAD2DEG*(theta-target_theta));
			cout << "angle: " << target_angle << " da: " << da << endl;
			cout << "l: " << l << " m: " << m << endl;
		}
	
		i = 2*m - l;
		j = m - 2*l;
		
		theta2 = atan((2*i-j)/(sqrt(3.0)*j));
		
		cout << "theta1: " << theta << " theta2: " << theta2 << " angle: " << RAD2DEG*theta << endl;
		
		Lx = abs(l*a*cos(theta) - 0.5*m*a*(sqrt(3.0)*sin(theta) + cos(theta)));
		Ly = abs(i*a*sin(theta) + 0.5*j*a*(sqrt(3.0)*cos(theta) - sin(theta)));
		Lz = dz_base;
		
		cout << "Lx: " << Lx << '\t' << "Ly: " << Ly << '\t' << "Lz: " << Lz <<endl;
			
		Nx = round(Lx/dx_base);
		Ny = round(Ly/dy_base);
		Nz = 1;
		
		cout << "Nx: " << Nx << '\t' << "Ny: " << Ny << '\t' << "Nz: " << Nz <<endl;
	}
		
	else if (IC_type == "HCP_ZROT"){
		
		l = 5; //sets the minimum size of the system
		
		cout << "Caluculating HCP 3D system size(Lx, Ly, Lz), rotated about z-axis" << endl; 
		
		while ((da > da_thresh) && (l<=50))
		{
			l = l+1;
			m = round(2.0*l*tan(target_theta)/(tan(target_theta)-sqrt(3.0)));
			theta = atan((sqrt(3.0)*m)/(m-2*l));
			da = abs(RAD2DEG*(theta-target_theta));
			cout << "angle: " << target_angle << " da: " << da << endl;
			cout << "l: " << l << " m: " << m << endl;
		}
	
		i = 2*m - l;
		j = m - 2*l;
		
		theta2 = atan((2*i-j)/sqrt(3.0));
		
		cout << "theta1: " << theta << " theta2: " << theta2 << " angle: " << RAD2DEG*theta << endl;
		
// 		alternative method for determining system size
// 		while (da > da_thresh)
// 		{
// 			l = l+1;
// 			k = round(0.5*l*(1-sqrt(3)*tan(target_theta)));
// 			theta = (atan((l-2*k)/(sqrt(3)*l)));
// 			da = abs(RAD2DEG*(theta-target_theta));
// 		}
// 	
// 		i = k - 2*l;
// 		j = 2*k - l;
// 		
// 		theta2 =  atan(sqrt(3.0)*j/(j-2*i));
		
		Lx = abs(l*a*cos(theta) - 0.5*m*a*(sqrt(3.0)*sin(theta) + cos(theta)));
		Ly = abs(i*a*sin(theta) + 0.5*j*a*(sqrt(3.0)*cos(theta) - sin(theta)));
		Lz = 4.0*c;
		
// 		Lx = 10.0*a;
// 		Ly = 10.0*sqrt(3.0)*a;
// 		Lz = 10.0*c;
		
		
		cout << "Lx: " << Lx << '\t' << "Ly: " << Ly << '\t' << "Lz: " << Lz <<endl;
	
		Nx = round(Lx/dx_base);
		Ny = round(Ly/dy_base);
		Nz = round(Lz/dy_base);
				
		cout << "Nx: " << Nx << '\t' << "Ny: " << Ny << '\t' << "Nz: " << Nz <<endl;
	}	
	
	else if (IC_type == "HCP_YROT"){
	
// 		double Lx1,Lz1;
// 		double Lx2,Lz2;
//		n3 = s_val; //sets the minimum size of the system
//		n1 = s_val; //sets the minimum size of the system
		cout << "Calculating HCP 3D system size(Lx, Ly, Lz), rotated about y-axis" << endl; 
		
		if (target_angle < 54.0) 
		{
		//cout << "if statement true" << endl;
			while ((da > da_thresh) && (n3<=200))
			{
				n3++;
				n1 = round(-(c/a)*n3*tan(target_theta));
				theta = atan(-(a/c)*(n1/n3));
				da = abs(RAD2DEG*(theta-target_theta));
				if (n3 == 200) {break;}
			}
		
			l1 = -8*n3;
			l3 = 3*n1;
			m1 = 1;
			theta2 = atan((c/a)*(l3/l1));

			cout << "n1: " << n1 << " n3: " << n3 << " l1: " << l1 << " l3: " << l3 << endl; 
			cout << "theta1: " << theta << " theta2: " << theta2 << " angle: " << RAD2DEG*theta << endl;
		
			Lx = abs(l1*a*cos(theta) + l3*c*sin(theta));
			Ly = sqrt(3)*a*2.0;
			Lz = abs(n3*c*cos(theta) - n1*a*sin(theta));
			
			cout << "Lx: " << Lx << '\t' << "Ly: " << Ly << '\t' << "Lz: " << Lz <<endl;
		
			Nx = round(Lx/dx_base);
			Ny = round(Ly/dy_base);
			Nz = round(Lz/dy_base);
		
			cout << "Nx: " << Nx << '\t' << "Ny: " << Ny << '\t' << "Nz: " << Nz <<endl;
		}	
		else
		{
			while ((da > da_thresh)&& (n3<=200))
			{
				n1++;
				n3 = round(-(a/c)*n1/(tan(target_theta)));
				theta = atan(-(a/c)*(n1/n3));
				da = abs(RAD2DEG*(theta-target_theta));
				if (n3 == 200) {break;
			}
		}
		
		l1 = -8*n3;
		l3 = 3*n1;
		m1 = 1;
		theta2 = atan((c/a)*(l3/l1));
		s_val = n3;

		cout << "n1: " << n1 << " n3: " << n3 << " l1: " << l1 << " l3: " << l3 << endl; 
		cout << "theta1: " << theta << " theta2: " << theta2 << " angle: " << RAD2DEG*theta << endl;
		
		Lx = abs(l1*a*cos(theta) + l3*c*sin(theta));
		Ly = sqrt(3)*a*2.0;
		Lz = abs(n3*c*cos(theta) - n1*a*sin(theta));
			
		cout << "Lx: " << Lx << '\t' << "Ly: " << Ly << '\t' << "Lz: " << Lz <<endl;
		
		Nx = round(Lx/dx_base);
		Ny = round(Ly/dy_base);
		Nz = round(Lz/dy_base);
		
		cout << "Nx: " << Nx << '\t' << "Ny: " << Ny << '\t' << "Nz: " << Nz <<endl;
		
		}
	}
	
	else if (IC_type == "HCP_UC"){

		l = 5; //sets the minimum size of the system
		
		//cout << "Caluculating HCP 3D system size(Lx, Ly, Lz), rotated about z-axis" << endl; 
		
		while ((da > da_thresh) && (l<=50))
		{
			l = l+1;
			m = round(2.0*l*tan(target_theta)/(tan(target_theta)-sqrt(3.0)));
			theta = atan((sqrt(3.0)*m)/(m-2*l));
			da = abs(RAD2DEG*(theta-target_theta));
			cout << "angle: " << target_angle << " da: " << da << endl;
			cout << "l: " << l << " m: " << m << endl;
		}
	
		i = 2*m - l;
		j = m - 2*l;
		
		theta2 = atan((2*i-j)/sqrt(3.0));
		
		//cout << "theta1: " << theta << " theta2: " << theta2 << " angle: " << RAD2DEG*theta << endl;
		
// 		alternative method for determining system size
// 		while (da > da_thresh)
// 		{
// 			l = l+1;
// 			k = round(0.5*l*(1-sqrt(3)*tan(target_theta)));
// 			theta = (atan((l-2*k)/(sqrt(3)*l)));
// 			da = abs(RAD2DEG*(theta-target_theta));
// 		}
// 	
// 		i = k - 2*l;
// 		j = 2*k - l;
// 		
// 		theta2 =  atan(sqrt(3.0)*j/(j-2*i));
		
		Lx = abs(l*a*cos(theta) - 0.5*m*a*(sqrt(3.0)*sin(theta) + cos(theta)));
		Ly = abs(i*a*sin(theta) + 0.5*j*a*(sqrt(3.0)*cos(theta) - sin(theta)));
		Lz = c;
		
		double UC_adjust = 1.0;

		if (target_angle == 0.0){
			Lx = Lx*UC_adjust;
			Ly = Ly*UC_adjust;
		}

		cout << "Lx: " << Lx << '\t' << "Ly: " << Ly << '\t' << "Lz: " << Lz <<endl;

		Nx = round(Lx/dx_base);
		Ny = round(Ly/dy_base);
		Nz = round(Lz/dy_base);

		cout << "Nx: " << Nx << '\t' << "Ny: " << Ny << '\t' << "Nz: " << Nz <<endl;

		// cout << "Setting HCP Unit cell size (no rotation)" << endl; 

		// c_num = 1.0;
		
		// double var = 0.00;
		
		// Lx = c_num*a;
		// Ly = c_num*sqrt(3.0)*a;
		// Lz = c_num*c;
			
		// Nx = round(Lx/dx_base);
		// Ny = round(Ly/dy_base);
		// Nz = round(Lz/dy_base);
		
		// Lx = Lx + var*Lx;
		// Ly = Ly + var*Ly;
		// Lz = Lz + var*Lz;		
		
// 		Lx = 5.0;
// 		Ly = 5.0*sqrt(3.0);
// 		Lz = 5.0*sqrt(8.0/3.0);

// 		//Zhen's parameters		
// 		Lx = 3.0;
// 		Ly = 3.0;
// 		Lz = 3.0;
// 		
// 		Nx = 32;
// 		Ny = 32;
// 		Nz = 32;
		
	}	
	
	else if (IC_type == "BCC_UC"){
		
		c_num = 3.0;
		
		a = sqrt(1.5);
			
		double var = 0.00;
		
		cout << "Setting BCC unit cell size (no rotation)" << endl;
		Lx = c_num*a;
		Ly = c_num*a;
		Lz = c_num*a;
		cout << "Lx: " << Lx << '\t' << "Ly: " << Ly << '\t' << "Lz: " << Lz <<endl;
	
		Nx = round(Lx/dx_base);
		Ny = round(Ly/dy_base);
		Nz = round(Lz/dy_base);
		
		if ((int)Nx % 2 == 1) {Nx = Nx + 1.0;}
		if ((int)Ny % 2 == 1) {Ny = Ny + 1.0;}
		if ((int)Nz % 2 == 1) {Nz = Nz + 1.0;}
		
// 		Nx = 128;
// 		Ny = 128;
// 		Nz = 128;

		Lx = Lx + var*Lx;
		Ly = Ly + var*Ly;
		Lz = Lz + var*Lz;
				
		cout << "Nx: " << Nx << '\t' << "Ny: " << Ny << '\t' << "Nz: " << Nz <<endl;
	}
	
	else if (IC_type == "FCC_UC"){
	
		cout << "Setting FCC unit cell size (no rotation)" << endl;
		
		c_num = 3.0;
		a = 1.5;
		
		double var = 0.007;
		
		Lx = c_num*a;
		Ly = c_num*a;
		Lz = c_num*a;
		
		cout << "Lx: " << Lx << '\t' << "Ly: " << Ly << '\t' << "Lz: " << Lz <<endl;
	
		Nx = round(Lx/dx_base);
		Ny = round(Ly/dy_base);
		Nz = round(Lz/dy_base);
		
		if ((int)Nx % 2 == 1) {Nx = Nx + 1.0;}
		if ((int)Ny % 2 == 1) {Ny = Ny + 1.0;}
		if ((int)Nz % 2 == 1) {Nz = Nz + 1.0;}
		
		Lx = Lx + var*Lx;
		Ly = Ly + var*Ly;
		Lz = Lz + var*Lz;
		
		cout << "Nx: " << Nx << '\t' << "Ny: " << Ny << '\t' << "Nz: " << Nz <<endl;
	}
	
	else if (IC_type == "SC_UC"){
	
		cout << "Setting FCC unit cell size (no rotation)" << endl;
		
		c_num = 2.0;
		
		Lx = c_num*a;
		Ly = c_num*a;
		Lz = c_num*a;
		
		cout << "Lx: " << Lx << '\t' << "Ly: " << Ly << '\t' << "Lz: " << Lz <<endl;
	
		Nx = round(Lx/dx_base);
		Ny = round(Ly/dy_base);
		Nz = round(Lz/dy_base);
		
		cout << "Nx: " << Nx << '\t' << "Ny: " << Ny << '\t' << "Nz: " << Nz <<endl;
	}
	
	else if (IC_type == "Columnar"){
	
		l = 20; //sets the minimum size of the system
		
		cout << "Caluculating HCP 3D system size for columnar system (z-axis rotation)" << endl; 
		
		while ((da > da_thresh) && (l<=50))
		{
			l = l+1;
			m = round(2.0*l*tan(target_theta)/(tan(target_theta)-sqrt(3.0)));
			theta = atan((sqrt(3.0)*m)/(m-2*l));
			da = abs(RAD2DEG*(theta-target_theta));
			cout << "angle: " << target_angle << " da: " << da << endl;
			cout << "l: " << l << " m: " << m << endl;
		}
	
		i = 2*m - l;
		j = m - 2*l;
		
		theta2 = atan((2*i-j)/sqrt(3.0));
		
		cout << "theta1: " << theta << " theta2: " << theta2 << " angle: " << RAD2DEG*theta << endl;
		
// 		alternative method for determining system size
// 		while (da > da_thresh)
// 		{
// 			l = l+1;
// 			k = round(0.5*l*(1-sqrt(3)*tan(target_theta)));
// 			theta = (atan((l-2*k)/(sqrt(3)*l)));
// 			da = abs(RAD2DEG*(theta-target_theta));
// 		}
// 	
// 		i = k - 2*l;
// 		j = 2*k - l;
// 		
// 		theta2 =  atan(sqrt(3.0)*j/(j-2*i));
		
		Lx = abs(l*a*cos(theta) - 0.5*m*a*(sqrt(3.0)*sin(theta) + cos(theta)));
		Ly = abs(i*a*sin(theta) + 0.5*j*a*(sqrt(3.0)*cos(theta) - sin(theta)));
		Lz = 2.0*c;
		
// 		Lx = 10.0*a;
// 		Ly = 10.0*sqrt(3.0)*a;
// 		Lz = 10.0*c;
		
		
		cout << "Lx: " << Lx << '\t' << "Ly: " << Ly << '\t' << "Lz: " << Lz <<endl;
	
		Nx = round(Lx/dx_base);
		Ny = round(Ly/dy_base);
		Nz = round(Lz/dy_base);
				
		cout << "Nx: " << Nx << '\t' << "Ny: " << Ny << '\t' << "Nz: " << Nz <<endl;
	}
	
		else if (IC_type == "Spherical"){
		double size_param;
		// l = 10; //sets the minimum size of the system
// 		
// 		cout << "Caluculating HCP 3D system size for columnar system (z-axis rotation)" << endl; 
// 		
// 		while ((da > da_thresh) && (l<=50))
// 		{
// 			l = l+1;
// 			m = round(2.0*l*tan(target_theta)/(tan(target_theta)-sqrt(3.0)));
// 			theta = atan((sqrt(3.0)*m)/(m-2*l));
// 			da = abs(RAD2DEG*(theta-target_theta));
// 			cout << "angle: " << target_angle << " da: " << da << endl;
// 			cout << "l: " << l << " m: " << m << endl;
// 		}
// 	
// 		i = 2*m - l;
// 		j = m - 2*l;
// 		
// 		theta2 = atan((2*i-j)/sqrt(3.0));
// 		
// 		cout << "theta1: " << theta << " theta2: " << theta2 << " angle: " << RAD2DEG*theta << endl;
		
// 		alternative method for determining system size
// 		while (da > da_thresh)
// 		{
// 			l = l+1;
// 			k = round(0.5*l*(1-sqrt(3)*tan(target_theta)));
// 			theta = (atan((l-2*k)/(sqrt(3)*l)));
// 			da = abs(RAD2DEG*(theta-target_theta));
// 		}
// 	
// 		i = k - 2*l;
// 		j = 2*k - l;
// 		
// 		theta2 =  atan(sqrt(3.0)*j/(j-2*i));
// 		
// 		Lx = abs(l*a*cos(theta) - 0.5*m*a*(sqrt(3.0)*sin(theta) + cos(theta)));
// 		Ly = abs(i*a*sin(theta) + 0.5*j*a*(sqrt(3.0)*cos(theta) - sin(theta)));
// 		Lz = 10.0*c;
		
        size_param = 20.0;

        Lx = size_param*a;
        Ly = size_param*sqrt(3.0)*a;
        Lz = size_param*c;
		
		
		cout << "Lx: " << Lx << '\t' << "Ly: " << Ly << '\t' << "Lz: " << Lz <<endl;
	
		Nx = round(Lx/dx_base);
		Ny = round(Ly/dy_base);
		Nz = round(Lz/dy_base);
				
		cout << "Nx: " << Nx << '\t' << "Ny: " << Ny << '\t' << "Nz: " << Nz <<endl;
	}
		
	else if (IC_type == "Liquid"){
	
		cout << "Setting FCC unit cell size (no rotation)" << endl;
		Lx = a;
		Ly = a;
		Lz = a;
		cout << "Lx: " << Lx << '\t' << "Ly: " << Ly << '\t' << "Lz: " << Lz <<endl;
	
		Nx = round(Lx/dx_base);
		Ny = round(Ly/dy_base);
		Nz = round(Lz/dy_base);
		
		cout << "Nx: " << Nx << '\t' << "Ny: " << Ny << '\t' << "Nz: " << Nz <<endl;
	}
		
 	dx = Lx/Nx;
 	dy = Ly/Ny;
 	dz = Lz/Nz;
	
	//dx = dx_base;
	//dy = dy_base;
	//dz = dz_base;
	
	cout << "dx: " << dx << '\t' << "dy: " << dy << '\t' << "dz: " << dz <<endl;
	cout << endl;
		
	if (IC_type == "HEXAGONAL") {dz = 0.5;}
	
	//cout << temp1 << temp2 << endl;
	
}
