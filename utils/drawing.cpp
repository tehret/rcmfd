// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.


#include "drawing.h"

void draw_line(double *igray, int a0, int b0, int a1, int b1, float value, int width, int height)
{

  int bdx,bdy;
  int sx,sy,dx,dy,x,y,z,l;

  bdx = width;
  bdy = height;

  if (a0 < 0) a0=0;
  else if (a0>=bdx) a0=bdx-1;

  if (a1<0)  a1=0;
  else  if (a1>=bdx)   a1=bdx-1;

  if (b0<0) b0=0;
  else if (b0>=bdy) b0=bdy-1;

  if (b1<0) 	b1=0;
  else if (b1>=bdy) b1=bdy-1;

  if (a0<a1) { sx = 1; dx = a1-a0; } else { sx = -1; dx = a0-a1; }
  if (b0<b1) { sy = 1; dy = b1-b0; } else { sy = -1; dy = b0-b1; }
  x=0; y=0;

  if (dx>=dy)
    {
      z = (-dx) / 2;
      while (abs(x) <= dx)
        {

          l =  (y+b0)*bdx+x+a0;

          igray[l] = value;

          x+=sx;
          z+=dy;
          if (z>0) { y+=sy; z-=dx; }

        }

    }
  else
    {
      z = (-dy) / 2;
      while (abs(y) <= dy) {

        l = (y+b0)*bdx+x+a0;
        igray[l] = value;

        y+=sy;
        z+=dx;
        if (z>0) { x+=sx; z-=dy; }
      }
    }

}


void draw_line(float *igray, int a0, int b0, int a1, int b1, float value, int width, int height)
{

  int bdx,bdy;
  int sx,sy,dx,dy,x,y,z,l;

  bdx = width;
  bdy = height;

  if (a0 < 0) a0=0;
  else if (a0>=bdx) a0=bdx-1;

  if (a1<0)  a1=0;
  else  if (a1>=bdx)   a1=bdx-1;

  if (b0<0) b0=0;
  else if (b0>=bdy) b0=bdy-1;

  if (b1<0) 	b1=0;
  else if (b1>=bdy) b1=bdy-1;

  if (a0<a1) { sx = 1; dx = a1-a0; } else { sx = -1; dx = a0-a1; }
  if (b0<b1) { sy = 1; dy = b1-b0; } else { sy = -1; dy = b0-b1; }
  x=0; y=0;

  if (dx>=dy)
    {
      z = (-dx) / 2;
      while (abs(x) <= dx)
        {

          l =  (y+b0)*bdx+x+a0;

          igray[l] = value;

          x+=sx;
          z+=dy;
          if (z>0) { y+=sy; z-=dx; }

        }

    }
  else
    {
      z = (-dy) / 2;
      while (abs(y) <= dy) {

        l = (y+b0)*bdx+x+a0;
        igray[l] = value;

        y+=sy;
        z+=dx;
        if (z>0) { x+=sx; z-=dy; }
      }
    }

}


void draw_square(float *igray, int a0, int b0, int w0, int h0, float value, int width, int height) //Mariano Rodríguez
{
        draw_line(igray,a0,b0,a0+w0,b0,value,width,height);
        draw_line(igray,a0,b0,a0,b0+h0,value,width,height);
        draw_line(igray,a0+w0,b0,a0+w0,b0+h0,value,width,height);
        draw_line(igray,a0,b0+h0,a0+w0,b0+h0,value,width,height);
}


void draw_parallelograms(float *igray, int* a, int* b, int* c, int* d, float value, int width, int height) //Mariano Rodríguez
{
        draw_line(igray,a[0],b[0],a[1],b[1],value,width,height);
        draw_line(igray,b[0],c[0],b[1],c[1],value,width,height);
        draw_line(igray,c[0],d[0],c[1],d[1],value,width,height);
        draw_line(igray,c[0],a[0],c[1],a[1],value,width,height);
}
