// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.


#ifndef UTILS_H
#define UTILS_H

#include <math.h>

void combine(float *u,float a,float *v,float b, float *w,  int size);
void compute_gradient_orientation(float* igray,float *grad, float *ori, int width, int height);
void sample(float *im, float *out, float factor, int width, int height, int channels);

#endif // UTILS_H
