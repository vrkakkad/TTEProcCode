template <class Type>
struct arg_struct {
    int nrhs, flag, frame, frameOffset;
    mwSize sz1, sz2, old_sz1, old_sz2, osz1, osz2;
    Type *in, *in1, *u, *v, *out, min_in, min_u, max_u, min_v, max_v;
};

template <class Type>
void *thread_execute(void *ptr) {
    arg_struct *p = (struct arg_struct *)ptr;
    int i, j, m;
    if (p.nrhs == 7) {
        // Search for the minimum value
        p.min_in = p.in[p.frame*p.old_sz1*p.old_sz2];
        for (i = 1; i < p.old_sz1*p.old_sz2; i++)
            p.min_in = min(p.in[i+p.frame*p.old_sz1*p.old_sz2], p.min_in);
    }
    
    // Pad image on all sides with the minimum value
    for(i = 0; i < p.sz1; i++) {
        for(j = 0; j < p.sz2; j++) {
            if( i == 0 || j == 0 || i == p.sz1-1 || j == p.sz2-1)
                p.in1[j*p.sz1+i] = p.min_in;
            else
                p.in1[j*p.sz1+i] = p.in[(j-1)*p.old_sz1+(i-1)+p.frame*p.old_sz1*p.old_sz2];
        }
    }
    
    // Use nearest neighbor method
    if(p.flag == 1) {
        for(i = 0; i < p.osz1*p.osz2; i++)
            p.v[i] = floor(p.v[i]);
        for(i = 0; i < p.osz1; i++) {
            for(j = 0; j < p.osz2; j++) {
                m = j*p.osz1+i;
                p.out[i*p.osz2+j + p.frameOffset] = p.in1[ (int)p.v[m]*sz1 + (int)p.u[m] ];
            }
        }
    }
    
    // Use (bi)linear interpolation
    // Since data is highly oversampled radially (grid spacing vs. actual
    // radial sampling intervals), can do nearest neighbor assumption for
    // radial portion. interpolate in theta direction
    else {
        Type err;
        for(i = 0; i < p.osz1; i++) {
            for(j = 0; j < p.osz2; j++) {
                m = j*p.osz1+i;
                err = p.v[m] - floor(p.v[m]);
                if(err != 0.0f)
                    p.out[i*p.osz2+j + p.frameOffset] = p.in1[ (int)floor(p.v[m])*sz1 + (int)p.u[m] ] +
                            p.in1[ (int) ceil(p.v[m])*sz1 + (int)p.u[m] ] -
                            p.in1[ (int)floor(p.v[m])*sz1 + (int)p.u[m] ])*err;
            }
        }
    }
    // Now just fix sector shape
    for(i = 0; i < p.osz1; i++) {
        for(j = 0; j < p.osz2; j++) {
            m = j*p.osz1+i;
            if(p.u[m] == p.min_u || p.u[m] == p.max_u || p.v[m] == p.min_v || p.v[m] == p.max_v)
                p.out[i*p.osz2+j + p.frameOffset] = p.min_in;
        }
    }
    
}

template <class Type>
void sector(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Declare input variables //
    Type *in;
    double min_phi, span_phi, apex, fs, *boxparams;
    int flag;
    
    // Declare other variables //
    double max_phi, incx, incy, x_min, x_max, y_min, y_max;
    int i, j, c, m, ndims, frame, frameOffset;
    const mwSize *in_dims;
    mwSize *out_dims, sz1, sz2, sz3, old_sz1, old_sz2, osz1, osz2;
    Type *in1, *x, *y, *r, *phi, *u, *v;
    Type *out, *ax, *lat;
    Type min_in, min_u, max_u, min_v, max_v;
        
    // Initialize input variables //
    in       = (Type *)mxGetData(prhs[1]);
    sz1      = mxGetM(prhs[1]);
    sz2      = mxGetN(prhs[1]);  // These are corrected later if necessary
    sz3      = 1;                // These are corrected later if necessary
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
        sz2 = in_dims[1];
        for(i = 2; i < ndims; i++)
            sz3 *= in_dims[i];
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
                mexErrMsgTxt("Argument #8, the value of points outside of the sector scan, must be a scalar.");
        }
    }
    else if (nrhs == 9) {
        if (mxGetM(prhs[7]) == 1 && mxGetN(prhs[7]) == 1)
            min_in = mxGetScalar(prhs[7]);
        else
            mexErrMsgTxt("Argument #8, the value of points outside of the sector scan, must be a scalar.");
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
    osz1 = floor((y_max-y_min)/incy)+1;
    osz2 = floor((x_max-x_min)/incx)+1;
    
    // Pre-allocate memory
    in1 = (Type *)mxMalloc(sizeof(Type)* sz1* sz2);
    x   = (Type *)mxMalloc(sizeof(Type)*osz1*osz2);
    y   = (Type *)mxMalloc(sizeof(Type)*osz1*osz2);
    r   = (Type *)mxMalloc(sizeof(Type)*osz1*osz2);
    phi = (Type *)mxMalloc(sizeof(Type)*osz1*osz2);
    u   = (Type *)mxMalloc(sizeof(Type)*osz1*osz2);
    v   = (Type *)mxMalloc(sizeof(Type)*osz1*osz2);
    
    // Determine dimensions of output
    out_dims = (mwSize *)mxMalloc(sizeof(mwSize)*ndims);
    out_dims[0] = osz2;
    out_dims[1] = osz1;
    for (i = 2; i < ndims; i++)
        out_dims[i] = in_dims[i];

    // Allocate output matrix
    plhs[0] = mxCreateNumericArray(ndims, out_dims, mxGetClassID(prhs[1]), mxREAL);
    plhs[1] = mxCreateNumericMatrix(1, osz2, mxGetClassID(prhs[1]), mxREAL);
    plhs[2] = mxCreateNumericMatrix(1, osz1, mxGetClassID(prhs[1]), mxREAL);
    
    // Set up output variables
    out = (Type *)mxGetData(plhs[0]);
    ax  = (Type *)mxGetData(plhs[1]);
    lat = (Type *)mxGetData(plhs[2]);

    // Calculate ax and lat and convert to mm
    for(i = 0; i < osz2; i++)
        ax[i] = (x_min+incx*i+(Type)apex)*100;
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

    // Set up multithreading
    numCPU = sysconf( _SC_NPROCESSORS_ONLN )
    pthread_t threads[numCPU];
    int thread_status[numCPU];
    
    struct arg_struct *p[numCPU];
    for (i = 0; i < numCPU; i++) {
        p[i].nrhs        = nrhs;
        p[i].flag        = flag;
        p[i].frame       = frame;
        p[i].frameOffset = frameOffset;
        p[i].sz1         = sz1;
        p[i].sz2         = sz2;
        p[i].old_sz1     = old_sz1;
        p[i].old_sz2     = old_sz2;
        p[i].osz1        = osz1;
        p[i].osz2        = osz2;
        p[i].min_in      = min_in;
        p[i].min_u       = min_u;
        p[i].max_u       = max_u;
        p[i].min_v       = min_v;
        p[i].max_v       = max_v;
        
        p[i].u   = mxMalloc(sizeof(Type)*osz1*osz2);
        p[i].v   = mxMalloc(sizeof(Type)*osz1*osz2);
        p[i].in  = mxMalloc(sizeof(Type)*old_sz1*old_sz2);
        p[i].in1 = mxMalloc(sizeof(Type)*sz1*sz2);
        p[i].out = mxMalloc(sizeof(Type)*osz1*osz2);
        
        for (j = 0; j < osz1*osz2; j++) {
            p[i].u[j] = u[j];
            p[i].v[j] = v[j];
        }
    }
    // Iterate through frames
    for (frame = 0; frame < sz3; frame++) {
        frameOffset = frame*osz1*osz2;
        
//         if (nrhs == 7) {
//             // Search for the minimum value
//             min_in = in[frame*old_sz1*old_sz2];
//             for (i = 1; i < old_sz1*old_sz2; i++)
//                 min_in = min(in[i+frame*old_sz1*old_sz2], min_in);
//         }
//         
//         // Pad image on all sides with the minimum value
//         for(i = 0; i < sz1; i++) {
//             for(j = 0; j < sz2; j++) {
//                 if( i == 0 || j == 0 || i == sz1-1 || j == sz2-1)
//                     in1[j*sz1+i] = min_in;
//                 else
//                     in1[j*sz1+i] = in[(j-1)*old_sz1+(i-1)+frame*old_sz1*old_sz2];
//             }
//         }
//         
//         // Use nearest neighbor method
//         if(flag == 1) {
//             for(i = 0; i < osz1*osz2; i++)
//                 v[i] = floor(v[i]);
//             for(i = 0; i < osz1; i++) {
//                 for(j = 0; j < osz2; j++) {
//                     m = j*osz1+i;
//                     out[i*osz2+j + frameOffset] = in1[ (int)v[m]*sz1 + (int)u[m] ];
//                 }
//             }
//         }
//         
//         // Use (bi)linear interpolation
//         // Since data is highly oversampled radially (grid spacing vs. actual
//         // radial sampling intervals), can do nearest neighbor assumption for
//         // radial portion. interpolate in theta direction
//         else {
//             Type err;
//             for(i = 0; i < osz1; i++) {
//                 for(j = 0; j < osz2; j++) {
//                     m = j*osz1+i;
//                     err = v[m] - floor(v[m]);
//                     if(err != 0.0f)
//                         out[i*osz2+j + frameOffset] = in1[ (int)floor(v[m])*sz1 + (int)u[m] ] +
//                                 in1[ (int) ceil(v[m])*sz1 + (int)u[m] ] -
//                                 in1[ (int)floor(v[m])*sz1 + (int)u[m] ])*err;
//                 }
//             }
//         }  
//         // Now just fix sector shape
//         for(i = 0; i < osz1; i++) {
//             for(j = 0; j < osz2; j++) {
//                 m = j*osz1+i;
//                 if(u[m] == min_u || u[m] == max_u || v[m] == min_v || v[m] == max_v)
//                     out[i*osz2+j + frameOffset] = min_in;
//             }
//         }
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