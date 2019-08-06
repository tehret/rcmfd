// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.

#include "flimage.h"



//////////////////////////////////////////////// Class flimage
//// Construction
flimage::flimage() : width(0), height(0), channels(0), p(0), pgray(0) 
{

}	 

flimage::flimage(int w, int h) : width(w), height(h), channels(1), p(new float[w*h]), pgray(p) 
{
	for (int j=width*height-1; j>=0 ; j--) p[j] = 0.0;
}	

flimage::flimage(int w, int h, int c) : width(w), height(h), channels(c), p(new float[w*h*c]), pgray(new float[w*h])
{
	for (int j=width*height*channels-1; j>=0 ; j--) p[j] = 0.0;
	for (int j=width*height-1; j>=0 ; j--) pgray[j] = 0.0;
}	

flimage::flimage(int w, int h, float v) : width(w), height(h), channels(1), p(new float[w*h]), pgray(p)
{
	for (int j=width*height-1; j>=0 ; j--) p[j] = v;
}

flimage::flimage(int w, int h, int c, float v) : width(w), height(h), channels(c), p(new float[w*h*c]), pgray(new float[w*h])
{
	for (int j=width*height*channels-1; j>=0 ; j--) p[j] = v;
	for (int j=width*height-1; j>=0 ; j--) pgray[j] = v;
}

flimage::flimage(int w, int h, float* v) : width(w), height(h), channels(1), p(new float[w*h]), pgray(p)
{
	for (int j=width*height-1; j>=0 ; j--) p[j] = v[j];
}

flimage::flimage(int w, int h, int c, float* v) : width(w), height(h), channels(c), p(new float[w*h*c]), pgray(new float[w*h])
{
    for(int i = 0; i < height; ++i) 
    for(int j = 0; j < width; ++j) 
    {
        float w = 0.;
        for(int c = 0; c < channels; ++c) 
        {
            p[i*width+j+c*width*height] = v[i*width+j+c*width*height];
            w += v[i*width+j+c*width*height];
        }
        pgray[i*width+j] = w/channels;
    }
}

void flimage::create(int w, int h, int c)
{ 
 	erase();
	width = w; height = h; channels = c;
	p = new float[w*h*c];	
	for (int j=width*height*channels-1; j>=0 ; j--) p[j] = 0.0;
	pgray = new float[w*h];	
	for (int j=width*height-1; j>=0 ; j--) pgray[j] = 0.0;
}

void flimage::create(int w, int h, int c, float* v)
{
 	erase();
	width = w; height = h; channels = c;  
    p = new float[w*h*c]; 
    pgray = new float[w*h]; 
    for(int i = 0; i < height; ++i) 
    for(int j = 0; j < width; ++j) 
    {
        float w = 0.;
        for(int c = 0; c < channels; ++c) 
        {
            p[i*width+j+c*width*height] = v[i*width+j+c*width*height];
            w += v[i*width+j+c*width*height];
        }
        pgray[i*width+j] = w/channels;
    }
}


flimage::flimage(const flimage& im) : width(im.width), height(im.height), channels(im.channels), p(new float[im.width*im.height*im.channels]), pgray(new float[im.width*im.height]) 
{
	for (int j=width*height*channels-1; j>=0 ; j--) p[j] = im.p[j];
	for (int j=width*height-1; j>=0 ; j--) pgray[j] = im.pgray[j];
}

flimage&  flimage::operator= (const flimage& im)
{	
	if (&im == this) {
		return *this;
	}
	
	if (width != im.width || height != im.height || channels != im.channels)
	{  			
	  	erase();
		width = im.width; height = im.height; channels = im.channels; 
        p = new float[width*height*channels];
        pgray = new float[width*height];
	}
	
	for (int j=width*height*channels-1; j>=0 ; j--) p[j] = im.p[j];
	for (int j=width*height-1; j>=0 ; j--) pgray[j] = im.pgray[j];
	return *this;	
}


//// Destruction
void flimage::erase() 
{
	width = height = channels = 0; 
    if(pgray != p) delete[] pgray;
    pgray = 0;
	if (p) delete[] p;
	p=0;
} 

void flimage::update(float* im) 
{

    for(int i = 0, k = 0; i < height; ++i) 
    for(int j = 0; j < width; ++j) 
    {
        float v = 0.;
        for(int c = 0; c < channels; ++c, ++k)
        {
            p[i*width+j+c*width*height] = im[i*width+j+c*width*height];
            v += im[i*width+j+c*width*height];
        }
        pgray[i*width+j] = v/channels;
    }
} 

flimage::~flimage()
{
	erase();
}




