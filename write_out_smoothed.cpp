#include "header.h"
#include <blitz/array.h>
#include <fstream>
#include <string.h>

extern "C"
{
	#include "visit_writer.h"
}

using namespace blitz;
using namespace std;

void GB_system::write_out_smooth(Array<double,3> array, string filename) 
{	
//  cout << "Real only VTK writer" << endl;
	//if(filename == "C"){
		ofstream myfile;
  		myfile.open(filename + ".txt");
  		myfile.precision(16);
  		myfile.setf(std::ios::fixed, std:: ios::floatfield); // floatfield set to fixed
		//float max_val = 0.0, min_val = 0.0;
  		//myfile << Ny << endl << Nx << endl << Nz << endl;
  		for (int i = 0; i<Nz; i++){
  			for (int j=0; j<Ny; j++){
  				for (int k=0; k<Nx; k++){
  				myfile << i << "\t" << j << "\t" << k << "\t" << array(i,j,k) << endl;
  				}
  			//myfile << endl;
  			}
  		//myfile << endl;	
  		}
	myfile.close();
	//}
    int NX = (int)Nx;
    int NY = (int)Ny;
    int NZ = (int)Nz;
    
    cout << "Writing " + filename << endl;  
    
    float temp = 0.0;
    
    int dims[] = {NX,NY,NZ};
    int nvars = 1;
    int vardims[] = {1};
    int centering[] = {1};
    const char *varnames[] = {"value"};
    float value[NZ][NY][NX];
    float *vars[] = { (float *)value};
    
    //cout << filename << ' ' << array.size() << endl;
    
    for (int i = 0; i < NZ; i++){
    	for (int j = 0; j < NY; j++){
    		for (int k = 0; k < NX; k++){
    			temp = array(i,j,k);
    			if((filename == "C") && (temp <= 10e-20)) {temp = 0.0;}
    			value[i][j][k] = temp;
    			//cout << temp << " " << value[i][j][k] << endl;
    		}
    	}
    }
    write_regular_mesh((filename + ".vtk").c_str(), 0, dims, nvars, vardims, centering, varnames, vars);
}	
//  cout << "Real only VTK writer" << endl;
	// if(filename == "C"){
// 		ofstream myfile;
//   		myfile.open(filename + ".txt");
//   		myfile.precision(16);
//   		myfile.setf(std::ios::fixed, std:: ios::floatfield); // floatfield set to fixed
// 		//float max_val = 0.0, min_val = 0.0;
//   		//myfile << Ny << endl << Nx << endl << Nz << endl;
//   		for (int i = 0; i<Nz; i++){
//   			for (int j=0; j<Ny; j++){
//   				for (int k=0; k<Nx; k++){
//   				myfile << i << "\t" << j << "\t" << k << "\t" << array(i,j,k) << endl;
//   				}
//   			//myfile << endl;
//   			}
//   		//myfile << endl;	
//   		}
// 	myfile.close();
// 	}
//     int NX = (int)Nx;
//     int NY = (int)Ny;
//     int NZ = (int)Nz;
//     
//     cout << "Writing " + filename << endl;  
//     
//     float temp_c = 0.0;
//     float temp_u = 0.0;
//     float temp_d = 0.0;
//     float temp_l = 0.0;
//     float temp_r = 0.0;
//     float temp_t = 0.0;
//     float temp_b = 0.0;
//     
//     int dims[] = {NX,NY,NZ};
//     int nvars = 1;
//     int vardims[] = {1};
//     int centering[] = {1};
//     const char *varnames[] = {"value"};
//     float value[NZ][NY][NX];
//     float *vars[] = { (float *)value};
//     
//     //cout << filename << ' ' << array.size() << endl;
//     
//     for (int i = 0; i < NZ; i++){
//     	for (int j = 0; j < NY; j++){
//     		for (int k = 0; k < NX; k++){
//     			temp_c = array(i,j,k);
//     			if ((i > 0) && (i < (NZ-1)) && (j > 0) && (j < (NY-1)) && (k > 0) && (k < (NX-1))){
//     				//temp_t = array(i,j,k);
//     				//temp_b = array(i,j,k-1);
//     				//temp_u = array(i,j+1,k);
//     				//temp_d = array(i,j-1,k);
//     				//temp_r = array(i+1,j,k);
//     				//temp_l = array(i-1,j,k);
//     				temp_c = (temp_c + temp_u + temp_d + temp_l + temp_r + temp_t + temp_b);
//     				}
//     			value[i][j][k] = temp_c;
//     			//cout << temp << " " << value[i][j][k] << endl;
//     		}
//     	}
//     }
//     write_regular_mesh((filename + ".vtk").c_str(), 0, dims, nvars, vardims, centering, varnames, vars);
// }

void GB_system::write_out_smooth(Array<complex<double>,3> array, string filename) 
{	
 cout << "Complex VTK writer disabled - no file written" << endl;
// 	ofstream myfile;
//   	myfile.open(filename);
//   	myfile.precision(16);
//   	myfile.setf(std::ios::fixed, std:: ios::floatfield); // floatfield set to fixed
//   	
//   	ofstream myfile_complex;
//   	myfile_complex.open("complex_" + filename);
//   	myfile_complex.precision(16);
//   	myfile_complex.setf(std::ios::fixed, std:: ios::floatfield); // floatfield set to fixed
//   	
//   	myfile << Ny << endl << Nx << endl << Nz << endl;
//   	myfile_complex << Ny << endl << Nx << endl << Nz << endl;
//   	for (int k = 0; k<Nz; k++){
//   		for (int j=0; j<Nx; j++){
//   			for (int i=0; i<Ny; i++){
//   				myfile << real(array(i,j,k)) << endl;
//   				myfile_complex << imag(array(i,j,k)) << endl;
//   			}
//   			//myfile<< endl;
//   			//myfile_complex << endl;
//   		}
//   		//myfile << endl;
//   		//myfile_complex << endl;
//   	}
//     myfile.close();
//     myfile_complex.close();
//     int NX = (int)Nx;
//     int NY = (int)Ny;
//     int NZ = (int)Nz;
//     float index = 0.0;
//     
//     cout << NX << NY << NZ << endl;
//     
//     int dims[] = {NX,NY,NZ};
//     int nvars = 3;
//     int vardims[] = {1,1,1};
//     int centering[] = {0,1};
//     const char *varnames[] = {"zonal", "nodal", "c_nodal"};
//     float zonal[NZ-1][NY-1][NX-1], nodal[NZ][NY][NX], c_nodal[NZ][NY][NX];
//     float *vars[] = {(float *)zonal, (float *)nodal, (float *)c_nodal};
//     for (int k = 0; k < NZ-1; ++k){
//     	for (int j = 0; j < NY-1; ++j){
//     		for (int i = 0; i < NX-1; ++i){
//     			++index;
//     			zonal[k][j][i] = index;
//     		}
//     	}
//     }
//     for (int k = 0; k < NZ-1; ++k){
//     	for (int j = 0; j < NY-1; ++j){
//     		for (int i = 0; i < NX-1; ++i){
//     			nodal[k][j][i] = real(array(j,i,k));
//     			c_nodal[k][j][i] = imag(array(j,i,k));
//     		}
//     	}
//     }    
//   	write_regular_mesh(filename.c_str(), 0, dims, nvars, vardims, centering, varnames, vars);
}


