#include "matlab_SMB_odefun.h"

void odefun(double *dydt, double t, double *y, int neq, const mxArray *args[])
{
    double tk;

    tk = mxGetScalar(args[0]);
    args = args+1;
    matlab_SMB_fun_y(tk,y,dydt);
}