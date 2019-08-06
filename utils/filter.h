// Credit: Rafael Grompone von Gioi


#ifndef FILTER_H
#define FILTER_H

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>

#define MAX(i,j) ( (i)<(j) ? (j):(i) )
#define MIN(i,j) ( (i)<(j) ? (i):(j) )

#define PI 3.14159265358979323846

float* gaussian_convolution(float *u, int width, int height, int channels, float sigma);

#endif // FILTER_H

