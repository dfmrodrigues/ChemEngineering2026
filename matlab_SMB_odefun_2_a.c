#include "matlab_SMB_odefun.h"

void odefun_2(double *dadt, double t, double *y, double *a, int neq, const mxArray *args[])
{
    double tk;
    int checkp;
    
    tk = mxGetScalar(args[0]);
    checkp = (mxGetScalar(args[1])!=0);
    args = args+2;
    matlab_SMB_fun_a(tk,a,y,dadt,checkp);
}