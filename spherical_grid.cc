#include "spherical_grid.h"
#include <iostream>
#include <cmath>


///////////////////////////////////////////////////////////////////////////////
// spherical_grid default constructor
///////////////////////////////////////////////////////////////////////////////
spherical_grid::spherical_grid(): m(pow(2,7)), R(1.){
    Ez = new float[ 2*m*m ];
    Hx = new float[ 2*m*(m-1) ];
    Hy = new float[ 2*m*m ];
    this->set_deltas();
    this->set_areas();
    this->set_fields_zero();
}
///////////////////////////////////////////////////////////////////////////////
// spherical simple constructor with net parameter
///////////////////////////////////////////////////////////////////////////////
spherical_grid::spherical_grid(int mm, float RR): m(mm), R(RR){
    Ez = new float[ 2*m*m ];
    Hx = new float[ 2*m*(m-1) ];
    Hy = new float[ 2*m*m ];
    this->set_deltas();
    this->set_areas();
    this->set_fields_zero();
}
///////////////////////////////////////////////////////////////////////////////
// default destructor
///////////////////////////////////////////////////////////////////////////////
spherical_grid::~spherical_grid(){
    delete [] Ez;
    delete [] Hx;
    delete [] Hy;
    delete [] area;
    delete [] delta_w_e;
}

///////////////////////////////////////////////////////////////////////////////
// functions computing sizes and areas of cells
///////////////////////////////////////////////////////////////////////////////
void spherical_grid::set_deltas(){
    delta_s_n = R * M_PI / m;
    delta_w_e = new float[m-1];
    for( int i = 0; i < m-1; ++i){
        delta_w_e[i] = sin( (i+1) * M_PI / m ) * R * M_PI / m;
    }

}

void spherical_grid::set_areas(){
    area = new float[m];
    area[0] = delta_w_e[0] * delta_s_n * sin(acos(delta_w_e[0]/(2.*delta_s_n))) /2.;
    for( int i = 1; i < m-1; ++i ){
        area[i] = delta_s_n * ( delta_w_e[i-1] + delta_w_e[i] )/2.;
    }
    area[m-1] = delta_w_e[m-2] * delta_s_n * sin(acos(delta_w_e[m-2]/(2.*delta_s_n))) /2.;
}

///////////////////////////////////////////////////////////////////////////////
// set initial zero values
///////////////////////////////////////////////////////////////////////////////
void spherical_grid::set_fields_zero(){
    for( int i = 0 ; i < 2*m*m; ++i ){
        Ez[i] = Hy[i] = 0.;
    }
   for( int i = 0 ; i < 2*m*(m-1); ++i ){
        Hx[i] = 0.;
    }
}

///////////////////////////////////////////////////////////////////////////////
// update Ez due to time step dt
///////////////////////////////////////////////////////////////////////////////
void spherical_grid::update_Ez( float dt, float epsilon0){
    for(int i = 0; i < 2*m-1; ++i){
        for(int j = 1; j < m-1; ++j){
            * get_Ez(i,j) += ( (*get_Hy(i+1,j) - *get_Hy(i,j))*delta_s_n 
                            + (*get_Hx(i,j-1))*delta_w_e[j-1] 
                            - (*get_Hx(i,j))*delta_w_e[j]
                               ) * dt/epsilon0/area[j];
        }
    }
    // glue side
    for(int j = 1; j < m-1; ++j){
        * get_Ez(2*m-1,j) += ( (*get_Hy(0,j) - *get_Hy(2*m-1,j))*delta_s_n 
                           + (*get_Hx(2*m-1,j-1))*delta_w_e[j-1] 
                           - (*get_Hx(2*m-1,j))*delta_w_e[j]
                              ) * dt/epsilon0/area[j];
        
    }
    // up glue
    for(int i = 0; i < 2*m-1; ++i){
        * get_Ez(i,m-1) += ( (*get_Hx(i,m-2))*delta_w_e[m-2] 
                            + (*get_Hy(i+1,m-1) - *get_Hy(i,m-1))*delta_s_n
                            )* dt/epsilon0/area[m-1];
    }
    // down glue
    for(int i = 0; i < 2*m-1; ++i){
        * get_Ez(i,0) += ( (*get_Hx(i,0))*delta_w_e[0] 
                            + (*get_Hy(i+1,0) - *get_Hy(i,0))*delta_s_n
                            )* dt/epsilon0/area[0];
    }
    
    // corners glue
    * get_Ez(2*m-1,m-1) += ( (*get_Hx(2*m-1,m-2))*delta_w_e[m-2] 
                            + (*get_Hy(0,m-1) - *get_Hy(2*m-1,m-1))*delta_s_n
                            )* dt/epsilon0/area[m-1];
    * get_Ez(2*m-1,0) += ( (*get_Hx(2*m-1,0))*delta_w_e[0] 
                            + (*get_Hy(0,0) - *get_Hy(2*m-1,0))*delta_s_n
                            )* dt/epsilon0/area[0];

}
///////////////////////////////////////////////////////////////////////////////
// update Hx due to time step dt
///////////////////////////////////////////////////////////////////////////////
void spherical_grid::update_Hx( float dt, float mu0){
    for(int i = 0; i < 2*m; ++i){
        for(int j = 0; j < m-1; ++j){
            *get_Hx(i,j) += ( (*get_Ez(i,j)) - (*get_Ez(i,j+1)) ) 
                            * dt/mu0/delta_w_e[j];
        }
    }
}
///////////////////////////////////////////////////////////////////////////////
// update Hy due to time step dt
///////////////////////////////////////////////////////////////////////////////
void spherical_grid::update_Hy( float dt, float mu0){
    for(int i = 1; i < 2*m; ++i){
        for(int j = 0; j < m; ++j){
            *get_Hy(i,j) += ((*get_Ez(i,j)) - (*get_Ez(i-1,j))) * dt/mu0/delta_s_n;
        }
    }
    for(int j = 0; j < m; ++j){
        *get_Hy(0,j) +=  ((*get_Ez(0,j))-(*get_Ez(2*m-1,j)))* dt/mu0/delta_s_n;
    }
}


///////////////////////////////////////////////////////////////////////////////
// prints
///////////////////////////////////////////////////////////////////////////////

void spherical_grid::snapshot_Ez(){
    std::cout << Ez[0] << " " << Hx[0] << " " << Hy[0] << std::endl;
    /*
    for(int i = 0; i < 2*m; ++i ){
        for(int j = 0; j < m; ++j){
           std::cout << i << " " << j << " " << *get_Hy(i,j) << std::endl; 
        }
    }
    std::cout << std::endl;
    std::cout << std::endl;
    */
}
void spherical_grid::print_areas(){
    for(int j = 0; j < m-1; ++j){
        std::cout << area[j] <<  " "  << delta_w_e[j] << std::endl;
    }
    std::cout << area[m-1] << " " << delta_s_n << std::endl;
}


