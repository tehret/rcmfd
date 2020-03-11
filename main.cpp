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


/**
 * @brief Creates an image showing the matches provided
 *
 * @param w,h,channels: size of the image (respectively width, height and channels)
 * @param im: image corresponding to the matches
 * @param matchings: list of matches to draw on the image
 **/
void write_images_matches(int channels, std::vector<float>& im,int w, int h, Matchingslist& matchings, string output)
{
    int sq = 2;

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
            draw_square(outimg[c],  round(matchings[i].first.x)-sq, round(matchings[i].first.y)-sq, 2*sq, 2*sq, value, w, h);
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
void write_image_mask(int ps, int w, int h, Matchingslist& matchings, string output)
{

    float * outmask = new float[w*h];

    for(int i=0; i < w*h; ++i)
        outmask[i] = 0;

    // Draw matches in red onto the grayscale output
    for(int i=0; i < (int) matchings.size(); i++)
    {

        // Draw the square of the first descriptor
        int sq = (int)(matchings[i].first.scale * ps/2);
        for(int xx = std::max((int)round(matchings[i].first.x)-sq, 0); xx < std::min((int)round(matchings[i].first.x)+sq, w); ++xx) 
        for(int yy = std::max((int)round(matchings[i].first.y)-sq, 0); yy < std::min((int)round(matchings[i].first.y)+sq, h); ++yy) 
            outmask[xx + w*yy] = 255.;

        // Draw the square of the second descripor
        sq = (int)(matchings[i].second.scale * ps/2);
        for(int xx = std::max((int)round(matchings[i].second.x)-sq, 0); xx < std::min((int)round(matchings[i].second.x)+sq, w); ++xx) 
        for(int yy = std::max((int)round(matchings[i].second.y)-sq, 0); yy < std::min((int)round(matchings[i].second.y)+sq, h); ++yy) 
            outmask[xx + w*yy] = 255.;
    }

    iio_save_image_float_vec(output.c_str(), outmask, w, h, 1);

    delete[] outmask;
}

//int maskConnectedComponent(unsigned entry, float th, float* mask, float* b1, float* b2, int w, int h, int c, std::vector<bool>& visited)
//{
//	std::stack<unsigned> stack;
//
//    int count = 1;
//    mask[entry] = 255;
//    visited[entry] = true;
//
//    stack.push(entry);
//    while(!stack.empty())
//    {
//        unsigned index = stack.top();
//        stack.pop();
//
//        unsigned a,b;
//        a = index % w;
//        b = index / w;
//
//        // Define the list of the neighbors used to compute the component
//        std::vector<std::pair<int,int>> neighbors{{a+1,b},{a-1,b},{a,b+1},{a,b-1}};
//
//        for(auto n : neighbors) 
//            // Add the possible neighbors to the current component
//            if((n.first >= 0) && (n.second >= 0) && (n.first < w) && (n.second < h))
//            {
//                bool test = true;
//                for(int ch = 0; ch < c; ++ch)
//                    if(std::abs(b1[n.first + w*n.second + ch*w*h] - b2[n.first + w*n.second + ch*w*h]) > th)
//                        test = false;
//                if(!visited[n.first + w*n.second] && test)
//                {
//                    mask[n.first + w*n.second] = 255;
//                    visited[n.first + w*n.second] = true;
//                    stack.push(n.first + w*n.second);
//                    count++;
//                }
//            }
//    }
//	return count;
//}
//
//void write_image_mask(int ps, double th, std::vector<float> image, int w, int h, int c, Matchingslist& matchings, string output)
//{
//    float * outmask = new float[w*h];
//    float * test = new float[w*h*c];
//
//    for(int i=0; i < w*h; ++i)
//        outmask[i] = 0;
//	std::vector<bool> visited(w*h, false);
//
//    printf("th %f\n", th);
//    for(int i=0; i < matchings.size(); ++i)
//    {
//        printf("doing %d out of %d\n", i, matchings.size());
//
//        // 2 -> 1
//        // blur im with scale1
//        float* blurred1 = image.data();
//        // blur im with scale2
//        float* blurred2 = gaussian_convolution(image.data(), w, h, c, matchings[i].second.scale);
//        // resample 2 on 1 
//        float* blurred2r = (float *) malloc( w * h * c * sizeof(float) );
//        float step =  matchings[i].second.scale / matchings[i].first.scale;
//        float theta = matchings[i].second.angle - matchings[i].first.angle;
//        float dx = step * cos(theta);
//        float dy = step * sin(theta);
//
//        for(int x = 0; x < w; ++x)
//        for(int y = 0; y < h; ++y)
//        {
//            float x_sample = dx * (x - matchings[i].first.x) - dy * (y - matchings[i].first.y) + matchings[i].second.x;
//            float y_sample = dy * (x - matchings[i].first.x) + dx * (y - matchings[i].first.y) + matchings[i].second.y;
//            for(int ch = 0; ch < c; ++ch)
//                blurred2r[x+y*w+ch*w*h] = interpolation(blurred2,w,h,x_sample,y_sample,ch);
//        }
//        // Compute the falsified region
//        maskConnectedComponent((int)(std::round(matchings[i].first.x) + w*std::round(matchings[i].first.y)), th, outmask, blurred1, blurred2r, w, h, c, visited);
//
//        
//        //free(blurred1);
//        free(blurred2);
//        free(blurred2r);
//
//        // 1 -> 2
//        // blur im with scale1
//        blurred1 = gaussian_convolution(image.data(), w, h, c, matchings[i].first.scale);
//        // blur im with scale2
//        blurred2 = gaussian_convolution(image.data(), w, h, c, matchings[i].second.scale);
//        // resample 2 on 1 
//        blurred2r = (float *) malloc( w * h * c * sizeof(float) );
//        step =  matchings[i].first.scale / matchings[i].second.scale;
//        theta = matchings[i].first.angle - matchings[i].second.angle;
//        dx = step * cos(theta);
//        dy = step * sin(theta);
//
//        for(int x = 0; x < w; ++x)
//        for(int y = 0; y < h; ++y)
//        {
//            float x_sample = dx * (x - matchings[i].second.x) - dy * (y - matchings[i].second.y) + matchings[i].first.x;
//            float y_sample = dy * (x - matchings[i].second.x) + dx * (y - matchings[i].second.y) + matchings[i].first.y;
//            for(int ch = 0; ch < c; ++ch)
//                blurred2r[x+y*w+ch*w*h] = interpolation(blurred2,w,h,x_sample,y_sample,ch);
//        }
//        // Compute the falsified region
//        maskConnectedComponent((int)(std::round(matchings[i].second.x) + w*std::round(matchings[i].second.y)), th, outmask, blurred1, blurred2r, w, h, c, visited);
//
//        
//        //free(blurred1);
//        //free(blurred2);
//        free(blurred2r);
//    }
//
//    iio_save_image_float_vec(output.c_str(), outmask, w, h, 1);
//    delete[] outmask;
//}

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
    double threshold = perform_matching(c, image, w, h, data, matchings, ps, tau, automatic);
    if (visual_output_path != "") 
        write_images_matches(c, image, w, h, matchings, visual_output_path);

    if (mask_output_path != "") 
        //write_image_mask(ps, threshold, image, w, h, c, matchings, mask_output_path);
        write_image_mask(ps, w, h, matchings, mask_output_path);

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
