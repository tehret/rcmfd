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

int main(int argc, char **argv)
{
	// Read parameters
	clo_usage("Copy-move forgery detection based on a contrario SIFT matching");
	clo_help(" NOTE: Input (<) and output (>) sequences are specified by their paths in printf format.\n");

	//! Paths to input/output sequences
	using std::string;
	string input_path = clo_option("-im"    , "" , "< Input image");
	string output_path = clo_option("-o"    , "data_matches.csv" , "> Output file containing the matches");
	string visual_output_path = clo_option("-vo"    , "output.png" , "> Output file containing the matches");

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
    perform_matching(c, image, w, h, data, matchings, ps, tau, automatic);
    write_images_matches(c, image, w, h, matchings, visual_output_path);

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
