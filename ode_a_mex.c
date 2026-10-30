#include "ode45_rp.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double t0, t, h, tnew, err, tfinal, htspan, threshold;
    int next, ntspan, haveeventfun;
    double rtol, absh, rh, hmax, hmin, tdir;
    double *tspan, temp;
    int done, nofailed, nfevals, nfailed, nsteps;
    int nout_new, oldnout, nout;
    double *tout, *aout;
    double *a0, *f0, *anew, *f, *a, *yi_int;
    double *v0, *v, *vnew, *direction, *indzc;
    int i, j, stop, neq;
    int nv;
    const mxArray **args;
    
    t0 = mxGetScalar(prhs[0]);
    a0 = mxGetPr(prhs[1]);
    neq = mxGetScalar(prhs[2]);
    tfinal = mxGetScalar(prhs[3]);
    htspan = mxGetScalar(prhs[4]);
    threshold = mxGetScalar(prhs[5]);
    nfevals = mxGetScalar(prhs[6]);
    nfailed = mxGetScalar(prhs[7]);
    rtol = mxGetScalar(prhs[8]);
    hmax = mxGetScalar(prhs[9]);
    tdir = mxGetScalar(prhs[10]);
    nsteps = mxGetScalar(prhs[11]);
    v0 = mxGetPr(prhs[12]);
    next = mxGetScalar(prhs[13]);
    ntspan = mxGetScalar(prhs[14]);
    haveeventfun = mxGetScalar(prhs[15]);
    tspan = mxGetPr(prhs[16]);
    nv = mxGetScalar(prhs[17]);
    yi_int = mxGetPr(prhs[18]);
    args = prhs+19;
    
    plhs[1] = mxCreateDoubleMatrix(neq,1,mxREAL);
    anew = mxGetPr(plhs[1]);
    plhs[3] = mxCreateDoubleMatrix(neq,7,mxREAL);
    f = mxGetPr(plhs[3]);
    plhs[11] = mxCreateDoubleMatrix(nv,1,mxREAL);
    vnew = mxGetPr(plhs[11]);
    plhs[12] = mxCreateDoubleMatrix(nv,1,mxREAL);
    direction = mxGetPr(plhs[12]);
    plhs[13] = mxCreateDoubleMatrix(nv,1,mxREAL);
    indzc = mxGetPr(plhs[13]);
    plhs[14] = mxCreateDoubleMatrix(nv,1,mxREAL);
    v = mxGetPr(plhs[14]);
    plhs[16] = mxCreateDoubleMatrix(1,ntspan,mxREAL);
    tout = mxGetPr(plhs[16]);
    plhs[17] = mxCreateDoubleMatrix(neq,ntspan,mxREAL);
    aout = mxGetPr(plhs[17]);
    plhs[20] = mxCreateDoubleMatrix(neq,1,mxREAL);
    a = mxGetPr(plhs[20]);
    
    t = t0;
    memcpy(a,a0,neq*sizeof(double));
    odefun_2(f+6*neq,t,yi_int+6*neq,a,neq,args);
    memcpy(vnew,v0,nv*sizeof(double));
    
    // hmin = 16*DBL_EPSILON*pow(2,floor(log(t0)*M_LOG2E));
    // absh = min(hmax,htspan);
    // rh = 0;
    // for(i = 0; i < neq; i++) {
    //     double rh_i;
    //     rh_i = f[i];
    //     rh_i = abs(rh_i)/max(abs(a[i]),threshold);
    //     if(rh_i > rh)
    //         rh = rh_i;
    // }
    // rh = rh/(0.8*pow(rtol,0.2));
    // if(absh*rh > 1)
    //     absh = 1/rh;
    // absh = max(absh,hmin);
    
    nout = 1;
    *tout = t;
    memcpy(aout,a,neq*sizeof(double));
    
    done = 0;
    while(!done) {
        
        // hmin = 16*DBL_EPSILON*pow(2,floor(log(t)*M_LOG2E));
        // absh = min(hmax, max(hmin, absh));
        // h = tdir * absh;
        // if(1.1*absh >= abs(tfinal - t)) {
        //     h = tfinal - t;
        //     absh = abs(h);
        //     done = 1;
        // }
        // else
        //     done = 0;
        if(tfinal == tspan[next-1])
            done = 1;
        else
            done = 0;
        nofailed = 1;
        next--;
        t = tspan[next];
        h = t-tspan[next-1];
        
        while(1) {
            tnew = t-h;
            for(i = 0; i < neq; i++) {
                anew[i] = 11.0/84.0*a[i];
            }
            odefun_2(f+5*neq,tnew,yi_int+(next-1)*7*neq+5*neq,anew,neq,args);
            tnew = t-h*8.0/9.0;
            for(i = 0; i < neq; i++) {
                anew[i] = -2187.0/6784.0*a[i];
                anew[i] = anew[i]-5103.0/18656.0*h*f[i+5*neq];
            }
            odefun_2(f+4*neq,tnew,yi_int+(next-1)*7*neq+4*neq,anew,neq,args);
            tnew = t-h*4.0/5.0;
            for(i = 0; i < neq; i++) {
                anew[i] = 125.0/192.0*a[i];
                anew[i] = anew[i]+49.0/176.0*h*f[i+5*neq];
                anew[i] = anew[i]-212.0/729.0*h*f[i+4*neq];
            }
            odefun_2(f+3*neq,tnew,yi_int+(next-1)*7*neq+3*neq,anew,neq,args);
            tnew = t-h*3.0/10.0;
            for(i = 0; i < neq; i++) {
                anew[i] = 500.0/1113.0*a[i];
                anew[i] = anew[i]+46732.0/5247.0*h*f[i+5*neq];
                anew[i] = anew[i]+64448.0/6561.0*h*f[i+4*neq];
                anew[i] = anew[i]+32.0/9.0*h*f[i+3*neq];
            }
            odefun_2(f+2*neq,tnew,yi_int+(next-1)*7*neq+2*neq,anew,neq,args);
            tnew = t-h*1.0/5.0;
            for(i = 0; i < neq; i++) {
                anew[i] = 0;
                anew[i] = anew[i]-355.0/33.0*h*f[i+5*neq];
                anew[i] = anew[i]-25360.0/2187.0*h*f[i+4*neq];
                anew[i] = anew[i]-56.0/15.0*h*f[i+3*neq];
                anew[i] = anew[i]+9.0/40.0*h*f[i+2*neq];
            }
            odefun_2(f+neq,tnew,yi_int+(next-1)*7*neq+neq,anew,neq,args);
            tnew = t;
            for(i = 0; i < neq; i++) {
                anew[i] = 35.0/384.0*a[i];
                anew[i] = anew[i]+9017.0/3168.0*h*f[i+5*neq];
                anew[i] = anew[i]+19372.0/6561.0*h*f[i+4*neq];
                anew[i] = anew[i]+44.0/45.0*h*f[i+3*neq];
                anew[i] = anew[i]+3.0/40.0*h*f[i+2*neq];
                anew[i] = anew[i]+1.0/5.0*h*f[i+1*neq];
            }
            odefun_2(f,tnew,yi_int+(next-1)*7*neq,anew,neq,args);
            for(i = 0; i < neq; i++) {
                anew[i] = a[i];
                anew[i] = anew[i]+h*f[i+5*neq];
                anew[i] = anew[i]+h*f[i+4*neq];
                anew[i] = anew[i]+h*f[i+3*neq];
                anew[i] = anew[i]+h*f[i+2*neq];
                anew[i] = anew[i]+h*f[i+1*neq];
                anew[i] = anew[i]+h*f[i+0*neq];
            }
            
            err = 0;
            
            nfevals = nfevals+6;
            
            if(err > rtol) {
                nfailed = nfailed + 1;
                if(nofailed) {
                    nofailed = 0;
                    absh = max(hmin,absh*max(0.1,0.8*pow(rtol/err,0.2)));
                }
                else
                    absh = max(hmin,0.5*absh);
                h = tdir * absh;
                done = 0;
            }
            else
                break;
        }
        nsteps++;
        stop = 0;
        
        if(haveeventfun) {
            memcpy(v,vnew,nv*sizeof(double));
            
            eventfun(tnew,anew,vnew,direction,nv,neq,args);
            
            stop = 0;
            for(i = 0; i < nv; i++) {
                indzc[i] = (sign(v[i])!=sign(vnew[i]))&&(direction[i]*(vnew[i]-v[i])>=0);
                if(indzc[i])
                    stop = 1;
            }
            if(stop)
                break;
        }
        
        nout_new = 0;
        nout_new++;
        next++;
        
        oldnout = nout;
        nout = nout+nout_new;
        next++;
        
        memcpy(tout+oldnout,&tnew,(nout-oldnout)*sizeof(double));
        memcpy(aout+oldnout*neq,anew,neq*(nout-oldnout)*sizeof(double));
        
        if(nofailed) {
            temp = 1.25*pow(err/rtol,0.2);
            if(temp > 0.2)
                absh = absh/temp;
            else
                absh = 5.0*absh;
        }
        t = tnew;
        memcpy(a,anew,neq*sizeof(double));
        memcpy(f+6*neq,f,neq*sizeof(double));
    }
    
    plhs[0] = mxCreateDoubleScalar(tnew);
    plhs[2] = mxCreateDoubleScalar(h);
    plhs[4] = mxCreateDoubleScalar(err);
    plhs[5] = mxCreateDoubleScalar(nfevals);
    plhs[6] = mxCreateDoubleScalar(nfailed);
    plhs[7] = mxCreateDoubleScalar(nofailed);
    plhs[8] = mxCreateDoubleScalar(absh);
    plhs[9] = mxCreateDoubleScalar(nsteps);
    plhs[10] = mxCreateDoubleScalar(stop);
    plhs[15] = mxCreateDoubleScalar(nout);
    plhs[18] = mxCreateDoubleScalar(next);
    plhs[19] = mxCreateDoubleScalar(t);
}