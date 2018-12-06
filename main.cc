#include<iostream>
#include<fstream>
#include<cmath>

#include "spherical_grid.h"


using namespace std;


int main(){
	
    ofstream ofs ("out.dat", ofstream::out);

    spherical_grid S(pow(2,7), 1.);

    float dt = 0.001;

    int t = 0;
    int t_max = 1000;

    
    while( t++ < t_max ){
        S.update_Ez(dt,1.);
        S.update_Hy(dt,1.);
        S.update_Hx(dt,1.);
        S.snapshot_Ez();
    }
    //S.print_areas();


	/*
	const double Re = 1.;
	const double epsilon0 = 1;
	const double mi0 = 1;
	
	const double dt = 0.001;
	const double time_steps = 1000;
	const int m = pow(2,8);
	
	const double d_theta = M_PI / m;
	const double d_fi = M_PI / m;
	const double delta_s_n = d_theta*Re;
	
	double delta_w_e[m-1];
	for(int j = 0; j < m-1; j++){
		delta_w_e[j] = Re * d_fi * sin((m-j-1)*M_PI/m);
	}	
	
	double S[m];
	S[m-1] = sin( acos( delta_w_e[m-2] / (2. * delta_s_n) ) ) * 
		    	delta_w_e[m-2] * delta_s_n / 2.;
	S[0] = sin( acos( delta_w_e[0] / (2. * delta_s_n) ) ) *
			delta_w_e[0] * delta_s_n / 2.;
	for(int j = 1; j < m-1; j++){
		S[j] = (delta_w_e[j-1] + delta_w_e[j]) * delta_s_n /2.;
	}	
	

    double ** Ez = new double * [2*m];
    for( int i=0; i< 2*m; ++i){
        Ez[i] = new double [m];
		for( int j = 0; j < m; j++ ){
            Ez[i][j] = 0.;
        }
    }
    double ** Hy = new double * [2*m];
    for( int i=0; i< 2*m; ++i){
        Hy[i] = new double [m];
		for( int j = 0; j < m; j++ ){
            Hy[i][j] = 0.;
        }
    }
    double ** Hx = new double * [2*m];
    for( int i=0; i< 2*m; ++i){
        Hx[i] = new double [m-1];
		for( int j = 0; j < m-1; j++ ){
            Hx[i][j] = 0.;
        }
    }

	//double Ez[2*m][m] = {};
	//double Hy[2*m][m] = {};
	//double Hx[2*m][m-1] = {};

	int t = 0;
	while( t++ < time_steps ){

///////////////////////////////////////////////////////////////////////////..
	// update Ez
    //
    if( t < 10 ){
    //Ez[0][0] = 1;
	}
	for( int i = 1; i < 2*m; i++ ){
		for( int j = 1; j < m-1; j++ ){
			Ez[i][j] = Ez[i][j] + (dt/(epsilon0*S[j]))* 
				( Hx[i][j-1] * delta_w_e[j-1] - Hx[i][j] * delta_w_e[j] + 
				delta_s_n * (Hy[i][j] - Hy[i-1][j]) );
		}
	}

	for( int i = 1; i < 2*m; i++ ){
		Ez[i][m-1] = Ez[i][m-1] + (dt/(epsilon0*S[m-1])) * 
				(  - Hx[i][m-2] * delta_w_e[m-2] + 
				delta_s_n * (Hy[i][m-1] - Hy[i-1][m-1]) );
		Ez[i][0] = Ez[i][0] + (dt/(epsilon0*S[0])) * 
				(  - Hx[i][0] * delta_w_e[0] + 
				delta_s_n * (Hy[i][0] - Hy[i-1][0]) );
			
	}

	for( int j = 1; j < m-1; j++ ){
		Ez[0][j] = Ez[0][j] + (dt/(epsilon0*S[j]))* 
			( Hx[0][j-1] * delta_w_e[j-1] - Hx[0][j] * delta_w_e[j] + 
			delta_s_n * (Hy[0][j] - Hy[2*m-1][j]) );
	}

	Ez[0][m-1] = Ez[0][m-1] + (dt/(epsilon0*S[m-1])) * 
				(  - Hx[0][m-2] * delta_w_e[m-2] + 
				delta_s_n * (Hy[0][m-1] - Hy[2*m-1][m-1]) );
	Ez[0][0] = Ez[0][0] + (dt/(epsilon0*S[0])) * 
			(  - Hx[0][0] * delta_w_e[0] + 
			delta_s_n * (Hy[0][0] - Hy[2*m-1][0]) );
	
	// update Hx
	
	for( int i = 0; i < 2*m; i++ ){
		for( int j = 0; j < m-1; j++){
			Hx[i][j] = Hx[i][j] + (dt/(mi0*delta_s_n)) *
					( Ez[i][j] - Ez[i][j+1] );

		}
	}
	
    // update Hy

	for( int i = 1; i < 2*m; i++ ){
		for( int j = 0; j < m; j++){
			Hy[i][j] = Hy[i][j] + (dt/(mi0*delta_w_e[j])) *
					( Ez[i][j] - Ez[i-1][j] );

		}
	}
	for( int j = 0; j < m; j++){
		Hy[0][j] = Hy[0][j] + (dt/(mi0*delta_w_e[j])) *
				( Ez[0][j] - Ez[2*m-1][j] );

	}
        ofs << t << " " << Ez[2*m-1][m-1] << endl;
    
///////////////////////////////////////////////////////////////////////////..
	}


    for(int i = 0; i < m; ++i) {
        delete [] Ez[i];
        delete [] Hy[i];
        delete [] Hx[i];
    }
    delete [] Ez;
    delete [] Hy;
    delete [] Hx;
*/    
    ofs.close();
	return 0;
}
