// Sector scan mode:
//
//       [out,ax,lat]=scan_convert('sector',in,min_phi,span_phi,apex,flag,fs,varargin);
// 
//       where   'sector' specifies sector scan mode,
//               in is the input data
//               min is the minimum (leftmost) sector angle (degrees)
//               span is the angle spaned by the sector (degrees)
//               apex is the distance in cm from the 'ducer face to the 
//                       center of curvature of the 'ducer.  Should be negative
//               flag specifies interpolation method (1 = nearest neighbor,
//                       2 = linear interpolation laterally)
//               fs specifies the axial sampling frequency
//               varargin{1} is an optional scalar argument to set the background value
//               varargin{2} is an optional vector argument for region selection,
//                       [ax_min, ax_max, ax_inc, lat_min, lat_max, lat_inc]
//
//
// Version 0.9
// 7/23/03  Stephen McAleavey, Duke University BME
//
// Revision History:
// 02/03/04 Brian Fahey, added flag variable to allow user to choose which geometric transformation to use
//       flag = 1 uses nearest neighbor
//       flag = 2 uses linear interpolation between theta angles
// radial oversampling makes nearest neighbor approx good for most pictures, although theta artifacts 
// can be seen if look closely. Interp method about 3x slower to implement. Recommend doing bulk processing
// with nearest neighbor, then making conference/paper images w/ interp
//
// 02/06/04 sam, added optional vector argument for region selection
//          [xmin xmax ymin ymax inc] all in meters
//          selects the rectangular window to zoom in on, as well as the grid 
// 	   spacing.  You can pick these numbers based on the axes of an 
//          un-zoomed image.  If this argument isn't supplied, the default
//          is the bounding rectangle of the sector. 
//
// 2011/05/04 Dongwoon Hyun
//      Mexified code.
#include "scan_convert.h"

void sectorSINGLE(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Declare input variables //
    float *in;
    double min_phi, span_phi, apex, fs, *boxparams;
    int flag;
    
    // Declare other variables //
    double max_phi, incx, incy, x_min, x_max, y_min, y_max;
    int i, j, c, m, ndims, frame, frameOffset;
    int sz1, sz2, sz3, old_sz1, old_sz2, osz1, osz2;
    const mwSize *in_dims;
    mwSize *out_dims;
    float *in1, *x, *y, *r, *phi, *u, *v;
    float *out, *ax, *lat;
    float min_in, min_u, max_u, min_v, max_v;
        
    // Initialize input variables //
    in       = (float *)mxGetPr(prhs[1]);
    sz1      = mxGetM(prhs[1]);
    sz2      = mxGetN(prhs[1]);
    sz3      = 1;
    in_dims  = mxGetDimensions(prhs[1]);
    ndims    = mxGetNumberOfDimensions(prhs[1]);
    min_phi  = mxGetScalar(prhs[2]);
    span_phi = mxGetScalar(prhs[3]);
    apex     = mxGetScalar(prhs[4]);
    flag     = (int)mxGetScalar(prhs[5]);
    fs       = mxGetScalar(prhs[6]);
    c = 1540;
    
    // Check to see how many frames there are
    if (ndims > 2) {
        sz2 = (int)in_dims[1];
        for(i = 2; i < ndims; i++)
            sz3 *= (int)in_dims[i];
    }
    
    // Increase sz1 and sz2 to pad image on all sides
    old_sz1 = sz1;
    old_sz2 = sz2;
    sz1 += 2;
    sz2 += 2;
    
    // Conversions to radians and meters //
    apex     = apex/100.0;
    min_phi  = min_phi *PI/180;
    span_phi = span_phi*PI/180;
    max_phi  = max(fabs(min_phi), fabs(min_phi + span_phi));
    
    // Use defaults for box parameters if not specified
    if (nrhs == 7 || nrhs == 8) {
        // Set scan conversion pitch //
        incx = 0.0002;
        incy = 0.0002;
        // Determine bounding box for scan data
        x_min = -cos(max_phi)*apex;
        x_max = sz1*(c/2)/fs-apex;
        y_min = min(0,x_max*sin(min_phi));
        y_max = max(0,x_max*sin(min_phi+span_phi));
        
        if (nrhs == 8) {
            if (mxGetM(prhs[7]) == 1 && mxGetN(prhs[7]) == 1)
                min_in = mxGetScalar(prhs[7]);
            else
                mexErrMsgTxt("Argument #8, min_in, must be a scalar.");
        }
    }
    else if (nrhs == 9) {
        if (mxGetNumberOfElements(prhs[8]) != 6)
            mexErrMsgTxt("Argument #9 must contain 6 elements: [x_min, x_max, incx, y_min, y_max, incy].");
        boxparams = mxGetPr(prhs[8]);
        x_min = boxparams[0];
        x_max = boxparams[1];
        incx  = boxparams[2];
        y_min = boxparams[3];
        y_max = boxparams[4];
        incy  = boxparams[5];
    }
    else
        mexErrMsgTxt("Incorrect number of input arguments");
    // Determine the size of the output
    osz1 = (int)floor((y_max-y_min)/incy)+1;
    osz2 = (int)floor((x_max-x_min)/incx)+1;
    
    // Pre-allocate memory
    in1 = (float *)mxMalloc(sizeof(float)* sz1* sz2);
    x   = (float *)mxMalloc(sizeof(float)*osz1*osz2);
    y   = (float *)mxMalloc(sizeof(float)*osz1*osz2);
    r   = (float *)mxMalloc(sizeof(float)*osz1*osz2);
    phi = (float *)mxMalloc(sizeof(float)*osz1*osz2);
    u   = (float *)mxMalloc(sizeof(float)*osz1*osz2);
    v   = (float *)mxMalloc(sizeof(float)*osz1*osz2);
    
    // Determine dimensions of output
    out_dims = (mwSize *)mxMalloc(sizeof(mwSize)*ndims);
    out_dims[0] = osz2;
    out_dims[1] = osz1;
    for (i = 2; i < ndims; i++)
        out_dims[i] = in_dims[i];

    // Allocate output matrix
    plhs[0] = mxCreateNumericArray(ndims, out_dims, mxSINGLE_CLASS, mxREAL);
    plhs[1] = mxCreateNumericMatrix(1, osz2, mxSINGLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericMatrix(1, osz1, mxSINGLE_CLASS, mxREAL);
    
    // Set up output variables
    out = (float *)mxGetPr(plhs[0]);
    ax  = (float *)mxGetPr(plhs[1]);
    lat = (float *)mxGetPr(plhs[2]);

    // Calculate ax and lat and convert to mm
    for(i = 0; i < osz2; i++)
        ax[i] = (x_min+incx*i+(float)apex)*100;
    for(i = 0; i < osz1; i++)
        lat[i]  = (y_min+incy*i)*100;
    
    // Create a grid of points in the output format
    // x,y is the output coordinate system
    // Equivalent to [x,y] = meshgrid(xmin:inc:x_max,y_min:inc:y_max);
    for(i = 0; i < osz1; i++) {
        for(j = 0; j < osz2; j++) {
            x[j*osz1+i] = x_min+j*incx;
            y[j*osz1+i] = y_min+i*incy;
        }
    }
    // Find these points in the transducer r-theta (polar) system
    for(i = 0; i < osz1*osz2; i++) {
        r[i]   = sqrt(x[i]*x[i]+y[i]*y[i]);
        phi[i] = atan2(y[i], x[i]);
    }
    // Turn distances into times to determine correct a-line sample
    // Turn angle into appropriate scan line
    for(i = 0; i < osz1*osz2; i++) {
        u[i] = (r[i]+apex)*fs/(c/2)-1;
        v[i] = ((phi[i]-min_phi)/span_phi)*(sz2-1);
    }
    // Limits to make sure you're in the input data
    // note that points outside the actual input data get mapped
    // to those border pixels created earlier
    for(i = 0; i < osz1*osz2; i++) {
        if(u[i] < 0.0f)
            u[i] = 0.0f;
        else if(u[i] > sz1-1)
            u[i] = sz1-1;
        if(v[i] < 0.0f)
            v[i] = 0.0f;
        else if(v[i] > sz2-1)
            v[i] = sz2-1;
        u[i] = floor(u[i]);
    }
    // Find max and min of u and v
    min_u = u[0];
    max_u = u[0];
    min_v = v[0];
    max_v = v[0];
    for(i = 1; i < osz1*osz2; i++) {
        if(u[i] > max_u)
            max_u = u[i];
        else if(u[i] < min_u)
            min_u = u[i];
        if(v[i] > max_v)
            max_v = v[i];
        else if(v[i] < min_v)
            min_v = v[i];
    }

    
    // Iterate through frames
    for (frame = 0; frame < sz3; frame++) {
        frameOffset = frame*osz1*osz2;
        if(frame%10 == 9) {
            mexPrintf("%d of %d frames completed...\n", frame+1, sz3);
            mexEvalString("drawnow"); 
        }
        if (nrhs == 7) {
            // Search for the minimum value
            min_in = in[frame*old_sz1*old_sz2];
            for (i = 1; i < old_sz1*old_sz2; i++)
                min_in = min(in[i+frame*old_sz1*old_sz2], min_in);
        }
        
        // Pad image on all sides with the minimum value
        for(i = 0; i < sz1; i++) {
            for(j = 0; j < sz2; j++) {
                if( i == 0 || j == 0 || i == sz1-1 || j == sz2-1)
                    in1[j*sz1+i] = min_in;
                else
                    in1[j*sz1+i] = in[(j-1)*old_sz1+(i-1)+frame*old_sz1*old_sz2];
            }
        }
        
        // Use nearest neighbor method
        if(flag == 1) {
            for(i = 0; i < osz1*osz2; i++)
                v[i] = floor(v[i]);
            #pragma omp parallel for num_threads(omp_get_num_procs()) collapse(2) shared(osz1, osz2, out, in1, v, sz1, u) private(j, m)
            for(i = 0; i < osz1; i++) {
                for(j = 0; j < osz2; j++) {
                    m = j*osz1+i;
                    out[i*osz2+j + frameOffset] = in1[ (int)v[m]*sz1 + (int)u[m] ];
                }
            }
        }
        
        // Use (bi)linear interpolation
        // Since data is highly oversampled radially (grid spacing vs. actual
        // radial sampling intervals), can do nearest neighbor assumption for
        // radial portion. interpolate in theta direction
        else {
            float err;
            #pragma omp parallel for num_threads(omp_get_num_procs()) collapse(2) shared(osz1, osz2, out, in1, v, sz1, u) private(j, m, err)
            for(i = 0; i < osz1; i++) {
                for(j = 0; j < osz2; j++) {
                    m = j*osz1+i;
                    err = v[m] - floor(v[m]);
                    if(err != 0.0f)
                        out[i*osz2+j + frameOffset] = in1[ (int)floor(v[m])*sz1 + (int)u[m] ] +
                                (in1[ (int)ceil(v[m])*sz1 + (int)u[m] ] -
                                in1[ (int)floor(v[m])*sz1 + (int)u[m] ])*err;
                }
            }
        }  
        // Now just fix sector shape
        for(i = 0; i < osz1; i++) {
            for(j = 0; j < osz2; j++) {
                m = j*osz1+i;
                if(u[m] == min_u || u[m] == max_u || v[m] == min_v || v[m] == max_v)
                    out[i*osz2+j + frameOffset] = min_in;
            }
        }
    }
        
    mxFree(u);
    mxFree(v);
    mxFree(r);
    mxFree(phi);
    mxFree(x);
    mxFree(y);
    mxFree(in1);
    mxFree(out_dims);
       
}
