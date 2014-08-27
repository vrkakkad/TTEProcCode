/*
 * ============================================================================
 * Name        : scan_convert.cpp
 * Author      : Dongwoon Hyun
 * ============================================================================
 *
 * ^ denotes option argument
 * For sector, input must be scan_convert(mode, in, min_phi, span_phi, apex, flag, fs, ^boxparams)
 * 
 */

#include "scan_convert.h"

// Gateway function for mex
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    // Argument check, determination of scan convert mode //
    if ( nrhs < 1 || !mxIsChar(prhs[0]) )
        mexErrMsgTxt("First argument must be the mode and must be a string.\n");
    char *mode;
    mwSize buflen;
    buflen = mxGetNumberOfElements(prhs[0]) + 1;
    mode   = (char *)mxMalloc(sizeof(char)*buflen);
    mxGetString(prhs[0], mode, buflen);
    
    // Convert all characters to lower case
    for(int i = 0; i != '\0'; i++) { mode[i] = tolower(mode[i]); }
    
    // Sector scan convert
    if (!strcmp(mode,"sector")) {
        // Argument Check //
        if (mxGetNumberOfElements(prhs[1]) == 0) { mexErrMsgTxt("Input matrix is empty."); }
        if (nrhs < 7)
            mexErrMsgTxt("Function is called by:\nscan_convert(mode,in,min_phi,span_phi,apex,flag,fs,extrapval,boxparams*);\n(* denotes optional argument)");
        if (!isScalar(prhs[2])) { mexErrMsgTxt("min_phi must be a scalar."); }
        if (!isScalar(prhs[3])) { mexErrMsgTxt("span_phi must be a scalar."); }
        if (!isScalar(prhs[4])) { mexErrMsgTxt("apex must be a scalar."); }
        if (!isScalar(prhs[5])) { mexErrMsgTxt("flag must be a scalar."); }
        if (!isScalar(prhs[6])) { mexErrMsgTxt("fs must be a scalar."); }
        if (nrhs >= 8 && !isScalar(prhs[7])) { mexErrMsgTxt("extrapval must be a scalar."); }
        if (nrhs > 9) { mexErrMsgTxt("Too many input arguments."); }
        
        // Call Function
        if (mxGetClassID(prhs[1]) == mxSINGLE_CLASS)
            sector<float >(nlhs, plhs, nrhs, prhs);
        else if (mxGetClassID(prhs[1]) == mxDOUBLE_CLASS)
            sector<double>(nlhs, plhs, nrhs, prhs);
        else
            mexErrMsgTxt("Input data must be of data type SINGLE or DOUBLE\n");
    }
//     else if (!strcmp(mode,"angle")) {
//         // Argument Check //
//         if (!( (nrhs == 2 && nlhs == 1) || (nrhs == 4 && nlhs == 3) ))
//             mexErrMsgTxt("Function is called by:\n[out,new_ax*,new_lat*]=scan_convert(mode,in,angle,axial*,lat*);\n(* denotes optional argument)");
//     }
}
// 
// void angle(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
//     // Declare input variables
//     float *in;
//     
// }

double max(double a, double b) { return a>b?a:b; }
float  max(float  a, float  b) { return a>b?a:b; }
double min(double a, double b) { return a<b?a:b; }
float  min(float  a, float  b) { return a<b?a:b; }
bool isScalar(mxArray *in) { return mxGetM(in) == 1 && mxGetN(in) == 1; }
bool isScalar(const mxArray *in) { return mxGetM(in) == 1 && mxGetN(in) == 1; }