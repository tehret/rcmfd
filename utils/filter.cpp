// Credit: Rafael Grompone von Gioi

#include "filter.h"

/*----------------------------------------------------------------------------*/
/* compute a Gaussian kernel of length n, standard deviation sigma,
   and centered at value mean.

   for example, if mean=0.5, the Gaussian will be centered in the middle point
   between values kernel[0] and kernel[1].

   kernel must be allocated to a size n.
 */
void gaussian_kernel(float * kernel, int n, float sigma, float mean)
{
  double sum = 0.0;
  double val;
  int i;

  /* check input */
  //if( kernel == NULL ) error("gaussian_kernel: kernel not allocated");
  //if( sigma <= 0.0 ) error("gaussian_kernel: sigma must be positive");

  /* compute Gaussian kernel */
  for(i=0; i<n; i++)
    {
      val = ( (float) i - mean ) / sigma;
      kernel[i] = exp( -0.5 * val * val );
      sum += kernel[i];
    }

  /* normalization */
  if( sum > 0.0 ) for(i=0; i<n; i++) kernel[i] /= sum;
}

/*----------------------------------------------------------------------------*/
/* filter an image with a Gaussian kernel of parameter sigma. return a pointer
   to a newly allocated filtered image, of the same size as the input image.
 */
float* gaussian_convolution(float * image, int X, int Y, int C, float sigma)
{
  int x,y,c,offset,i,j,nx2,ny2,n;
  float * kernel;
  float * tmp;
  float * out;
  float val,prec;

  /* check input */
  //if( sigma <= 0.0 ) error("gaussian_filter: sigma must be positive");
  //if( image == NULL || X < 1 || Y < 1 ) error("gaussian_filter: invalid image");

  /* get memory */
  out = (float *) malloc( X * Y * C * sizeof(float) );
  tmp = (float *) malloc( X * Y * sizeof(float) );

  /* compute gaussian kernel */
  /*
     The size of the kernel is selected to guarantee that the first discarded
     term is at least 10^prec times smaller than the central value. For that,
     the half size of the kernel must be larger than x, with
       e^(-x^2/2sigma^2) = 1/10^prec
     Then,
       x = sigma * sqrt( 2 * prec * ln(10) )
   */
  prec = 3.0;
  offset = (int) ceil( sigma * sqrt( 2.0 * prec * log(10.0) ) );
  n = 1 + 2 * offset; /* kernel size */
  kernel = (float *) malloc( n * sizeof(float) );
  gaussian_kernel(kernel, n, sigma, (float) offset);

  /* auxiliary variables for the double of the image size */
  nx2 = 2*X;
  ny2 = 2*Y;

  for(c=0; c<C; c++)
  {
      /* x axis convolution */
      for(x=0; x<X; x++)
          for(y=0; y<Y; y++)
          {
              val = 0.0;
              for(i=0; i<n; i++)
              {
                  j = x - offset + i;

                  /* symmetry boundary condition */
                  while(j<0) j += nx2;
                  while(j>=nx2) j -= nx2;
                  if( j >= X ) j = nx2-1-j;

                  val += image[j+y*X+c*X*Y] * kernel[i];
              }
              tmp[x+y*X] = val;
          }

      /* y axis convolution */
      for(x=0; x<X; x++)
          for(y=0; y<Y; y++)
          {
              val = 0.0;
              for(i=0; i<n; i++)
              {
                  j = y - offset + i;

                  /* symmetry boundary condition */
                  while(j<0) j += ny2;
                  while(j>=ny2) j -= ny2;
                  if( j >= Y ) j = ny2-1-j;

                  val += tmp[x+j*X] * kernel[i];
              }
              out[x+y*X+c*X*Y] = val;
          }
  }

  /* free memory */
  free( (void *) kernel );
  free( (void *) tmp );

  return out;
}
