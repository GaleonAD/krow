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

    void set_deltas();
    void set_areas();
    void set_fields_zero();

    public:
    // constructors
    spherical_grid();
    spherical_grid(int, float);
    ~spherical_grid();

    float * get_Ez(int i, int j){ return Ez + i*m+j; }
    float * get_Hx(int i, int j){ return Hx + i*(m-1)+j; }
    float * get_Hy(int i, int j){ return Hy + i*m+j; }
    void update_Ez(float, float);
    void update_Hx(float, float);
    void update_Hy(float, float);
    void snapshot_Ez();
    void print_areas();
    //spherical_grid(int M, float Re);
   

};


#endif
