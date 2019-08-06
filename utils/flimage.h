// Authors: Unknown. Please, if you are the author of this file, or if you 
// know who are the authors of this file, let us know, so we can give the 
// adequate credits and/or get the adequate authorizations.


#ifndef _FLIMAGE_H_
#define _FLIMAGE_H_

#include <iostream>
#include <string>
#include <cstdio>

class flimage {
	
private:
	
	int	width, height, channels;	// image size
	float*	p;		// array of color levels: level of pixel (x,y) is p[c*width*height+y*width+x]
	float*	pgray;	// array of gray levels: level of pixel (x,y) is p[y*width+x]
	
public:
	
	//// Construction
	flimage();
	flimage(int w, int h);
	flimage(int w, int h, float v);
	flimage(int w, int h, float* v);
	flimage(int w, int h, int c);
	flimage(int w, int h, int c, float v);
	flimage(int w, int h, int c, float* v);
	flimage(const flimage& im);
	flimage& operator= (const flimage& im);
	
	
	void create(int w, int h, int c);
	void create(int w, int h, int c, float *v);
	void create(int w, int h) {create(w,h,1);};
	void create(int w, int h, float *v) {create(w,h,1,v);};
	
	//// Destruction
	void erase();
	void update(float* im);
	~flimage();
	
	//// Get Basic Data	
	int nwidth() const {return width;} 	// image size
	int nheight() const {return height;} 
	int nchannels() const {return channels;} 
	
	/// Access values
    float* getColorPlane() const {return p;}	// return the adress of the array of values
    float* getPlane() const {return pgray;}	// return the adress of the array of values
	
	float operator()(int x, int y) const {return pgray[ y*width + x ];} 	// acces to the (x,y) value
	float& operator()(int x, int y) {return pgray[ y*width + x ];}	// by value (for const images) and by reference
	float operator()(int x, int y, int c) const {return p[c*width*height + y*width + x];} 	// acces to the (x,y) value
	float& operator()(int x, int y, int c) {return p[c*width*height + y*width + x ];}	// by value (for const images) and by reference
	
	
};


#endif

