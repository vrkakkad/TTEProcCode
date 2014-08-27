#include <stdio.h>
#include <stdlib.h>
#include "mex.h"
#include <math.h>
#include <string.h>
#include <iostream>

// Define constants
#define PI 3.14159265

// Prototypes of local functions
double max(double a, double b);
float  max(float  a, float  b);
double min(double a, double b);
float  min(float  a, float  b);
bool isScalar(mxArray *in);
bool isScalar(const mxArray *in);
// void angle (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

#include "sector.h"