// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.


#ifndef DRAWING_H
#define DRAWING_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>


#define MAX(i,j) ( (i)<(j) ? (j):(i) )
#define MIN(i,j) ( (i)<(j) ? (i):(j) )

// #define PI 3.14159
#define PI 3.14159265358979323846



void draw_line(float *igray, int a0, int b0, int a1, int b1, float value, int width, int height);

/**
 * @brief Draws a rectangle starting from the point \f$(a0,b0) \in [1,width]\times[1,height] \f$. The three other corners of the rectangle are: \f$(a0+w0,b0), (a0,b0+h0), (a0+w0,b0+h0).\f$
 * @param igray A one channel image onto which a rectangle is to be created.
 * @param (a0,b0) Image coordinates of the starting point of the rectangle
 * @param w0 width of the rectangle
 * @param h0 height of the rectangle
 * @param value Intensity of the lines belonging to the rectangle
 * @param width Image width
 * @param height Image height
 * @author Mariano Rodríguez
 */
void draw_square(float *igray, int a0, int b0, int w0, int h0, float value, int width, int height);

#endif
