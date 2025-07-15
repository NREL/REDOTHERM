#include "mex.h"
#include <math.h>

void panlener_dsfun(const double *delta, double *ds, mwSize n) {
    for (mwSize i = 0; i < n; ++i) {
        double d = delta[i] > 0 ? delta[i] : 1e-300;
        double x = log10(d);
        double val;

        if (x <= -3.6) {
            val = 291.1;
        } else if (x <= -1.03) {
            double x2 = x * x;
            double x3 = x2 * x;
            double x4 = x3 * x;
            double x5 = x4 * x;
            val = 40.1831 - 60.8925*x + 150.565*x2 + 126.019*x3 + 39.6964*x4 + 4.47333*x5;
        } else if (x <= -0.56) {
            double x2 = x * x;
            double x3 = x2 * x;
            val = -73.0801 - 726.348*x - 859.650*x2 - 368.734*x3;
        } else {
            val = 129.6;
        }

        ds[i] = val;
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
    double *ds = mxGetPr(plhs[0]);

    panlener_dsfun(delta, ds, n);
}
