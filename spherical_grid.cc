#include "spherical_grid.h"
#include <cmath>

spherical_grid::spherical_grid(): m(pow(2,7)), R(1.){
    Ez = new float[ 2*m*m ];
    Hx = new float[ 2*m*(m-1) ];
    Hy = new float[ 2*m*m ];
    this->compute_deltas();
    this->compute_areas();
    this->set_fields_zero();
}

void spherical_grid::compute_deltas(){
    delta_s_n = R * M_PI / m;
    delta_w_e = new float[m-1];
    for( int i = 0; i < m-1; ++i){
        delta_w_e[i] = sin( (i+1) * M_PI / m ) * R * M_PI / m;
    }

}

void spherical_grid::compute_areas(){
    area = new float[m];
    area[0] = delta_w_e[0] * delta_s_n * sin(acos(delta_w_e[0]/(2.*delta_s_n))) /2.;
    for( int i = 1; i < m-1; ++i ){
        area[i] = delta_s_n * ( delta_w_e[i-1] + delta_w_e[i] )/2.;
    }
    area[m-1] = delta_w_e[m-2] * delta_s_n * sin(acos(delta_w_e[m-2]/(2.*delta_s_n))) /2.;
}

void spherical_grid::set_fields_zero(){
    for( int i = 0 ; i < 2*m*m; ++i ){
        Ez[i] = Hy[i] = 0.;
    }
   for( int i = 0 ; i < 2*m*(m-1); ++i ){
        Hx[i] = 0.;
    }
}
