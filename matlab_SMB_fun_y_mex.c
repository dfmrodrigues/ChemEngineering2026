#include "matlab_SMB_odefun.h"
#include "ode45_rp.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double tk;
    double *y, *dydt;
    int M;

    tk = mxGetScalar(prhs[0]);
    y = mxGetPr(prhs[1]);
    M = mxGetM(prhs[1]);
    plhs[0] = mxCreateDoubleMatrix(M,1,mxREAL);
    dydt = mxGetPr(plhs[0]);
    matlab_SMB_fun_y(tk,y,dydt);
}