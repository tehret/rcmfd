// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.

#include "utils.h"

void combine(float *u,float a,float *v,float b, float *w,  int size)  { for(int i=0;i<size ;i++)   w[i]= a*u[i] + b*v[i];  };

void compute_gradient_orientation(float* igray,float *grad, float *ori, int width, int height)
{
    float xgrad, ygrad;
    int rows, cols, r, c;

    rows = height;
    cols = width;
    

    for (r = 0; r < rows; r++)
      for (c = 0; c < cols; c++) {
        if (c == 0)
          xgrad = 2.0 * (igray[r*cols+c+1] - igray[r*cols+c]);
        else if (c == cols-1)
          xgrad = 2.0 * (igray[r*cols+c] - igray[r*cols+c-1]);
        else
          xgrad = igray[r*cols+c+1] - igray[r*cols+c-1];
        if (r == 0)
          ygrad = 2.0 * (igray[r*cols+c] - igray[(r+1)*cols+c]);
        else if (r == rows-1)
          ygrad = 2.0 * (igray[(r-1)*cols+c] - igray[r*cols+c]);
        else
          ygrad = igray[(r-1)*cols+c] - igray[(r+1)*cols+c];
        

        if (grad) grad[r*cols+c] = (float)sqrt((double)(xgrad * xgrad + ygrad * ygrad));
        if (ori) ori[r*cols+c] = (float)atan2 (-(double)ygrad,(double)xgrad);
      
      }
};

void sample(float *im, float *out, float factor, int width, int height, int channels)
{

	int swidth = (int)((float) width / factor);
	int sheight = (int)((float) height / factor);
	
    for(int c=0; c < channels; ++c)
	for(int j=0; j < sheight; j++)
	 for(int i=0; i < swidth; i++)
		out[j*swidth + i + c*swidth*sheight] = im[(int)((float) j * factor) * width + (int) ((float) i*factor) + c*width*height];
};
