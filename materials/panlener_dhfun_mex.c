#include "mex.h"
#include <math.h>

void panlener_dhfun(const double *delta, double *dh, mwSize n) {
    for (mwSize i = 0; i < n; ++i) {
        double d = delta[i] > 0 ? delta[i] : 1e-300;
        double x = log10(d);
        double val;

        if (x <= -3.6) {
            val = 428.5;
        } else if (x <= -1.03) {
            double x2 = x * x;
            double x3 = x2 * x;
            double x4 = x3 * x;
            double x5 = x4 * x;
            val = 302.618 + 27.8543*x + 253.263*x2 + 184.585*x3 + 54.0575*x4 + 5.83517*x5;
        } else if (x <= -0.56) {
            double x2 = x * x;
            double x3 = x2 * x;
            double x4 = x3 * x;
            val = -643.250 - 5627.64*x - 10806.6*x2 - 8774.26*x3 - 2558.64*x4;
        } else {
            val = 408.5;
        }

        dh[i] = val*1e3;
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    if (nrhs != 1)
        mexErrMsgTxt("One input required.");
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
        mexErrMsgTxt("Input must be a real double array.");

    mwSize n = mxGetNumberOfElements(prhs[0]);
    const double *delta = mxGetPr(prhs[0]);

    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[0]), mxGetN(prhs[0]), mxREAL);
    double *dh = mxGetPr(plhs[0]);

    panlener_dhfun(delta, dh, n);
}