#include "matlab_SMB_odefun.h"

void eventfun(double t, double *y, double *v, double *direction, int nv, int neq, const mxArray **args)
{
    int i;
    
    for(i = 0; i < nv; i++) {
        v[i] = 0;
        direction[i] = 1.0;
    }
}
