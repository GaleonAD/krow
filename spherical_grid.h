#ifndef SPHERICAL_GRID_H
#define SPHERICAL_GRID_H


class spherical_grid{
    private:
    int m;
    float R;
    float * Ez;
    float * Hx;
    float * Hy;

    float delta_s_n;
    float * delta_w_e;
    float * area;

    void compute_deltas();
    void compute_areas();
    void set_fields_zero();

    public:
    spherical_grid();
    //spherical_grid(int M, float Re);
   

};


#endif
