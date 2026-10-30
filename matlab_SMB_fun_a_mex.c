#include "matlab_SMB_odefun.h"
#include "ode45_rp.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double tk;
    double *a, *y, *dadt;
    int M;

    tk = mxGetScalar(prhs[0]);
    a = mxGetPr(prhs[1]);
    M = mxGetM(prhs[1]);
    y = mxGetPr(prhs[2]);
    plhs[0] = mxCreateDoubleMatrix(M,1,mxREAL);
    dadt = mxGetPr(plhs[0]);
    matlab_SMB_fun_a(tk,a,y,dadt,1);
}