#ifndef CELL_H
#define CELL_H
#include "wall.h"

class cell{
    private:
    double Ez;

    wall * eleft;
    wall * eright;
    wall * eup;
    wall * edown;

    public:
    cell();

};

#endif
