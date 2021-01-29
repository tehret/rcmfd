/*
 * Copyright (c) 2019, Thibaud Ehret <ehret.thibaud@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU Affero General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <vector>
#include <stack>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include "matching.h"
#include "utils/cmd_option.h"
#include "utils/drawing.h"
extern "C"{
#include "utils/iio.h"
}
#include "utils/filter.h"
#include "utils/sift.h"
#include "inverse_compositional_algorithm.h"

#define PSIZE 49
#define PSIZE2 49

/**
 * @brief Creates an image showing the matches provided
 *
 * @param w,h,channels: size of the image (respectively width, height and channels)
 * @param im: image corresponding to the matches
 * @param matchings: list of matches to draw on the image
 **/
void write_images_matches(int channels, std::vector<float>& im, int ps, int w, int h, Matchingslist& matchings, string output)
{
    int sq;

    std::vector<float *> outimg;
    for(int c=0;c<3;c++)
        outimg.push_back(new float[w*h]);

    // Transform the image to grayscale and copy it to output
    for(int c=0;c<3;c++)
    for(int j = 0; j < (int) h; j++)
    for(int i = 0; i < (int) w; i++)
        outimg[c][j*w+i] = (channels > 1) ? (im[j*w+i] + im[1*w*h+j*w+i] + im[2*w*h+j*w+i])/3. : im[j*w+i];

    // Draw matches in red onto the grayscale output
    float value;
    for(int i=0; i < (int) matchings.size(); i++)
    for(int c=0;c<3;c++)
    {
            value =  (c == 0) ? 255. : 0.;
            // Draw a line between the two matches
            draw_line(outimg[c],  round(matchings[i].first.x), round(matchings[i].first.y),
                      round(matchings[i].second.x), round(matchings[i].second.y), value, w, h);

            // Draw both squares corresponding the position of the keypoints (not taking into account the scale nor the orientation)
            sq = (int)(matchings[i].first.scale * ps/2);
            draw_square(outimg[c],  round(matchings[i].first.x)-sq, round(matchings[i].first.y)-sq, 2*sq, 2*sq, value, w, h);
            sq = (int)(matchings[i].second.scale * ps/2);
            draw_square(outimg[c],  round(matchings[i].second.x)-sq, round(matchings[i].second.y)-sq, 2*sq, 2*sq, value, w, h);
    }

    // Reorganize data before saving the image
    float * rgb = new float[w*h*3];
    for(int c=0;c<3;c++)
    for(int j = 0; j < (int) h; j++)
    for(int i = 0; i < (int) w; i++)
        rgb[j*w*3+i*3+c] = outimg[c][j*w+i];
    iio_save_image_float_vec(output.c_str(), rgb, w, h, 3);

    delete[] rgb;
    for(int c=0;c<3;c++)
        delete[] outimg[c];
}

/**
 * @brief Creates a mask image (outdated very simple regions)
 *
 * @param w,h: size of the image (respectively width, height and channels)
 * @param matchings: list of matches to draw on the image
 **/
//void write_image_mask(int ps, int w, int h, Matchingslist& matchings, string output)
//{
//
//    float * outmask = new float[w*h];
//
//    for(int i=0; i < w*h; ++i)
//        outmask[i] = 0;
//
//    // Draw matches in red onto the grayscale output
//    for(int i=0; i < (int) matchings.size(); i++)
//    {
//
//        // Draw the square of the first descriptor
//        int sq = (int)(matchings[i].first.scale * ps/2);
//        for(int xx = std::max((int)round(matchings[i].first.x)-sq, 0); xx < std::min((int)round(matchings[i].first.x)+sq, w); ++xx) 
//        for(int yy = std::max((int)round(matchings[i].first.y)-sq, 0); yy < std::min((int)round(matchings[i].first.y)+sq, h); ++yy) 
//            outmask[xx + w*yy] = 255.;
//
//        // Draw the square of the second descripor
//        sq = (int)(matchings[i].second.scale * ps/2);
//        for(int xx = std::max((int)round(matchings[i].second.x)-sq, 0); xx < std::min((int)round(matchings[i].second.x)+sq, w); ++xx) 
//        for(int yy = std::max((int)round(matchings[i].second.y)-sq, 0); yy < std::min((int)round(matchings[i].second.y)+sq, h); ++yy) 
//            outmask[xx + w*yy] = 255.;
//    }
//
//    iio_save_image_float_vec(output.c_str(), outmask, w, h, 1);
//
//    delete[] outmask;
//}

int maskConnectedComponent(unsigned entry, float th, float* mask, float* b1, float* b2, int w, int h, int c, int* visited, int ref)
{
	std::stack<unsigned> stack;

    float varim = 0;
    float vardiff = 0;

    // Compute the variance of the patches
    unsigned a,b;
    a = entry % w;
    b = entry / w;
    for(int x = a-PSIZE2/2; x <= a+PSIZE2/2; ++x)
    for(int y = b-PSIZE2/2; y <= b+PSIZE2/2; ++y)
    for(int ch = 0; ch < c; ++ch)
    {
        varim += b1[x + y*w + ch*w*h];
        vardiff += b2[x + y*w + ch*w*h];
    }
    varim /= PSIZE2*PSIZE2*c;
    vardiff /= PSIZE2*PSIZE2*c;

    if(vardiff < varim)
    {
        stack.push(entry);
        for(int x = a-PSIZE2/2; x <= a+PSIZE2/2; ++x)
        for(int y = b-PSIZE2/2; y <= b+PSIZE2/2; ++y)
            mask[x + y*w] = 255;
    } 

    while(!stack.empty())
    {
        unsigned index = stack.top();
        stack.pop();

        a = index % w;
        b = index / w;

        // Define the list of the neighbors used to compute the component
        std::vector<std::pair<int,int>> neighbors{{a+PSIZE2/2+1,b},{a-PSIZE2/2-1,b},{a,b+PSIZE2/2+1},{a,b-PSIZE2/2-1}};

        for(auto n : neighbors) 
        {
            if((n.first-PSIZE2/2 >= 0) && (n.second-PSIZE2/2 >= 0) && (n.first+PSIZE2/2 < w) && (n.second+PSIZE2/2 < h))
                if(visited[n.first + w*n.second] <= ref)
                {
                    visited[n.first + w*n.second] = ref+1;
                    if(mask[n.first + w*n.second] == 0)
                    {
                        
                        varim = 0;
                        vardiff = 0;

                        for(int x = a-PSIZE2/2; x <= a+PSIZE2/2; ++x)
                            for(int y = b-PSIZE2/2; y <= b+PSIZE2/2; ++y)
                                for(int ch = 0; ch < c; ++ch)
                                {
                                    varim += b1[x + y*w + ch*w*h];
                                    vardiff += b2[x + y*w + ch*w*h];
                                }
                        varim /= PSIZE2*PSIZE2*c;
                        vardiff /= PSIZE2*PSIZE2*c;

                        if(vardiff < varim)
                        {
                            stack.push(entry);
                            for(int x = n.first-PSIZE2/2; x <= n.first+PSIZE2/2; ++x)
                                for(int y = n.second-PSIZE2/2; y <= n.second+PSIZE2/2; ++y)
                                    mask[x + y*w] = 255;
                            stack.push(n.first + w*n.second);
                        } 
                    }
                    else
                        stack.push(n.first + w*n.second);

                }
        }
    }
}

void write_image_mask(int ps, double th, std::vector<float> image, int w, int h, int c, Matchingslist& matchings, string output)
{
    float * outmask = new float[w*h];
    int * visited = new int[w*h];

    for(int i=0; i < w*h; ++i)
    {
        outmask[i] = 0;
        visited[i] = 0;
    }

    float* imgrad = (float *) malloc(w * h * c * sizeof(float) );
    float* resamp = (float *) malloc(w * h * c * sizeof(float) );
    float* diffgrad = (float *) malloc(w * h * c * sizeof(float) );
    float ddiff;
    for(int x = 1; x < w-1; ++x)
    for(int y = 1; y < h-1; ++y)
    for(int ch = 0; ch < c; ++ch)
            imgrad[x+y*w+ch*w*h] = (ddiff= (0.5 * (image[(x+1)+y*w+ch*w*h] - image[(x-1)+y*w+ch*w*h])))*ddiff + (ddiff= (0.5 * (image[x+(y+1)*w+ch*w*h] - image[x+(y-1)*w+ch*w*h])))*ddiff;

    for(int i=0; i < matchings.size(); ++i)
    {
        printf("Mask for match %d\n", i);
        // Extract patches
        float step =  matchings[i].second.scale / matchings[i].first.scale;
        double* patch1 = (double *) malloc(PSIZE * PSIZE * c * sizeof(double) );
        double* patch2 = (double *) malloc(PSIZE * PSIZE * c * sizeof(double) );
        bool flipped = matchings[i].second.flipped;
        for(int x = 0; x < PSIZE; ++x)
        for(int y = 0; y < PSIZE; ++y)
        for(int ch = 0; ch < c; ++ch)
        {
                patch1[ch+x*c+y*PSIZE*c] = image[(x-PSIZE/2 + (int)std::round(matchings[i].first.x))+(y-PSIZE/2 + (int)std::round(matchings[i].first.y))*w+ch*w*h] - matchings[i].first.mean[ch];
                patch2[ch+x*c+y*PSIZE*c] = interpolation(image.data(),w,h, (flipped?(PSIZE/2 - x-1):(x-PSIZE/2))*step + (int)std::round(matchings[i].second.x),(y-PSIZE/2)*step + (int)std::round(matchings[i].second.y),ch) - matchings[i].second.mean[ch];
        }

        // Register the local patches
        int nparams = 3;
        double* p = new double[3];
        p[0]=0;p[1]=0;p[2]=0;
        pyramidal_inverse_compositional_algorithm(
           patch1, patch2, p, nparams, PSIZE, PSIZE, c,
           2, 0.5, 0.001, 3, 0., 0, 1,
           5, 3, 0
        );

        // resample 2 on 1 
        float theta = p[2];
        float dx = step * cos(theta);
        float dy = step * sin(theta);

        for(int x = 0; x < w; ++x)
        for(int y = 0; y < h; ++y)
        {
            float x_sample;
            if(flipped)
                x_sample = dx * (matchings[i].first.x + PSIZE/2 - x-1) - dy * (y - matchings[i].first.y + PSIZE/2) + matchings[i].second.x - step*(PSIZE/2) - p[0];
            else
                x_sample = dx * (x-matchings[i].first.x + PSIZE/2) - dy * (y - matchings[i].first.y + PSIZE/2) + matchings[i].second.x - step*(PSIZE/2) + p[0];
            float y_sample = dy * (flipped?(matchings[i].first.x + PSIZE/2 - x-1):(x-matchings[i].first.x + PSIZE/2)) + dx * (y - matchings[i].first.y + PSIZE/2) + matchings[i].second.y - step*(PSIZE/2) + p[1];
            for(int ch = 0; ch < c; ++ch)
                resamp[x+y*w+ch*w*h] = interpolation(image.data(),w,h,x_sample,y_sample,ch);
        }

        for(int x = 1; x < w-1; ++x)
        for(int y = 1; y < h-1; ++y)
        for(int ch = 0; ch < c; ++ch)
                diffgrad[x+y*w+ch*w*h] = (ddiff= (0.5 * (image[(x+1)+y*w+ch*w*h] - image[(x-1)+y*w+ch*w*h]) - 0.5 * (resamp[(x+1)+y*w+ch*w*h] - resamp[(x-1)+y*w+ch*w*h])))*ddiff + (ddiff= (0.5 * (image[x+(y+1)*w+ch*w*h] - image[x+(y-1)*w+ch*w*h]) - 0.5 * (resamp[x+(y+1)*w+ch*w*h] - resamp[x+(y-1)*w+ch*w*h])))*ddiff;

        // Compute the falsified region
        maskConnectedComponent((int)(std::round(matchings[i].first.x) + w*std::round(matchings[i].first.y)), th, outmask, imgrad, diffgrad, w, h, c, visited, i);

        step =  matchings[i].first.scale / matchings[i].second.scale;
        for(int x = 0; x < PSIZE; ++x)
        for(int y = 0; y < PSIZE; ++y)
        for(int ch = 0; ch < c; ++ch)
        {
                patch1[ch+x*c+y*PSIZE*c] = image[(x-PSIZE/2 + (int)std::round(matchings[i].second.x))+(y-PSIZE/2 + (int)std::round(matchings[i].second.y))*w+ch*w*h] - matchings[i].second.mean[ch];
                patch2[ch+x*c+y*PSIZE*c] = interpolation(image.data(),w,h, (flipped?(PSIZE/2 - x-1):(x-PSIZE/2))*step + (int)std::round(matchings[i].first.x),(y-PSIZE/2)*step + (int)std::round(matchings[i].first.y),ch) - matchings[i].first.mean[ch];
        }

        // Register the local patches
        p[0]=0;p[1]=0;p[2]=0;
        pyramidal_inverse_compositional_algorithm(
           patch1, patch2, p, nparams, PSIZE, PSIZE, c,
           2, 0.5, 0.001, 3, 0., 0, 1,
           5, 3, 0
        );

        // resample 2 on 1 
        theta = p[2];
        dx = step * cos(theta);
        dy = step * sin(theta);

        for(int x = 0; x < w; ++x)
        for(int y = 0; y < h; ++y)
        {
            float x_sample;
            if(flipped)
                x_sample = dx * (matchings[i].second.x + PSIZE/2 - x-1) - dy * (y - matchings[i].second.y + PSIZE/2) + matchings[i].first.x - step*(PSIZE/2) - p[0];
            else
                x_sample = dx * (x-matchings[i].second.x + PSIZE/2) - dy * (y - matchings[i].second.y + PSIZE/2) + matchings[i].first.x - step*(PSIZE/2) + p[0];
            float y_sample = dy * (flipped?(matchings[i].second.x + PSIZE/2 - x-1):(x-matchings[i].second.x + PSIZE/2)) + dx * (y - matchings[i].second.y + PSIZE/2) + matchings[i].first.y - step*(PSIZE/2) + p[1];
            for(int ch = 0; ch < c; ++ch)
                resamp[x+y*w+ch*w*h] = interpolation(image.data(),w,h,x_sample,y_sample,ch);
        }

        for(int x = 1; x < w-1; ++x)
        for(int y = 1; y < h-1; ++y)
        for(int ch = 0; ch < c; ++ch)
                diffgrad[x+y*w+ch*w*h] = (ddiff= (0.5 * (image[(x+1)+y*w+ch*w*h] - image[(x-1)+y*w+ch*w*h]) - 0.5 * (resamp[(x+1)+y*w+ch*w*h] - resamp[(x-1)+y*w+ch*w*h])))*ddiff + (ddiff= (0.5 * (image[x+(y+1)*w+ch*w*h] - image[x+(y-1)*w+ch*w*h]) - 0.5 * (resamp[x+(y+1)*w+ch*w*h] - resamp[x+(y-1)*w+ch*w*h])))*ddiff;

        // Compute the falsified region
        maskConnectedComponent((int)(std::round(matchings[i].second.x) + w*std::round(matchings[i].second.y)), th, outmask, imgrad, diffgrad, w, h, c, visited, i);

    }

    iio_save_image_float_vec(output.c_str(), outmask, w, h, 1);
    delete[] outmask;
    free(imgrad);
    free(diffgrad);
    free(resamp);
}

int main(int argc, char **argv)
{
	// Read parameters
	clo_usage("Copy-move forgery detection based on a contrario SIFT matching");
	clo_help(" NOTE: Input (<) and output (>) sequences are specified by their paths in printf format.\n");

	//! Paths to input/output sequences
	using std::string;
	string input_path = clo_option("-im"    , "" , "< Input image");
	string output_path = clo_option("-o"    , "data_matches.csv" , "> Output file containing the matches");
	string visual_output_path = clo_option("-vo"    , "" , "> Output file containing the matches");
	string mask_output_path = clo_option("-mo"    , "" , "> Output file containing the mask");

	// Parameters for the matching
	int ps = clo_option("-ps", 8, "< Patch size for the descriptor");
	float tau = clo_option("-tau", 1., "< Threshold for the patch matching, also correspond to the NFA threshold when the automatic thershold is activated");
	bool grayscale = clo_option("-gs", false, "< Force the image to have only one channel");
	float automatic = clo_option("-auto"    , false, "< Use the threshold defined by NFA (using tau as parameter), otherwise use the provided threshold");

	if(input_path == "")
	{
		fprintf(stderr, "Options '-h', '-help' and '--help' list the different options available.\n");
		return EXIT_FAILURE;
	}

    // Load image, transform it into a one channel version if necessary
    std::vector<float> image;
    int w,h,c;
	float *imTmp;
	imTmp =  iio_read_image_float_vec(input_path.c_str(), &w, &h, &c);

    if(grayscale)
    {
        float *imTmpGS = new float[w*h];
        // Transform the input image into a single channel image
        for(int i = 0; i < w; ++i)
        for(int j = 0; j < h; ++j)
        {
            float temp = 0.;
            for(int k = 0; k < c; ++k)
                temp += imTmp[k + i*c + j*c*w];
            imTmpGS[i + j*w] = temp/(float)c;
        }
        free(imTmp);
        imTmp = imTmpGS;
        c = 1;
    }
    image.resize(w*h*c);
    
    for(int i = 0; i < w; ++i)
    for(int j = 0; j < h; ++j)
    for(int k = 0; k < c; ++k)
        image[i + j*w + k*w*h] = imTmp[k + i*c + j*c*w];
	free(imTmp);

    // Perform matching
    Matchingslist matchings;
    vector< float > data;
    double threshold = perform_matching(c, image, w, h, data, matchings, ps, tau, automatic, true);
    if (visual_output_path != "") 
        write_images_matches(c, image, ps, w, h, matchings, visual_output_path);

    if (mask_output_path != "") 
        write_image_mask(ps, threshold, image, w, h, c, matchings, mask_output_path);
        //write_image_mask(ps, w, h, matchings, mask_output_path);

    // Save results
    ofstream myfile;
    myfile.open(output_path.c_str(), std::ofstream::out | std::ofstream::trunc);
    myfile<<"x1, y1, sigma1, angle1, x2, y2, sigma2, angle2"<<endl;

    int wo = 8;
    if (matchings.size()>0)
    {
        int cont =1;
        myfile << ((double) data[0]) << ",";

        for ( int i = 1; i < (int) (wo*matchings.size()); i++ )
        {
            if (cont ==(wo-1))
            {
                myfile << ((double) data[i]) << endl;
                cont = 0;
            }
            else
            {
                myfile << ((double) data[i]) << ",";
                cont = cont +1;
            }

        }
    }
    myfile.close();

    return 0;
}
