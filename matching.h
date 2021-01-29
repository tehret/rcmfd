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

#ifndef MATCHING_H
#define MATCHING_H

#include <math.h>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <cstdio>
#include "utils/sift.h"
#include "stats/stats.hpp"

#define FILTER_RADIUS 4
#define NOTDEF -1024.0

/************* DATATYPES ********************/

/**
 * @brief Definition of a keypoint
 */
struct KeyPoint
{
    float x;
    float y;
    float scale;
    float angle;
    void* kp_ptr;
};

/**
 * @brief Used to aggregate spatially keypoints
 */
struct KeyPoints
{
    float x, y, sum_x, sum_y;
    std::vector<KeyPoint> KPvec;
};

/**
 * @brief Keypointlist is used to store keypoints during the computation
 */
typedef std::vector<KeyPoint> Keypointlist;

/**
 * @brief A simple version of a keypoint (without descriptor).
 */
struct Keypoint_simple {
    float x, y, scale, angle;
    float mean[3];
    bool flipped;
};

/**
 * @brief The very definition of a match as pair of keypoints.
 */
typedef std::pair<Keypoint_simple, Keypoint_simple> Matching;

/**
 * @brief List of matches.
 */
typedef std::vector<Matching> Matchingslist;

/********************************************/
/************* FUNCTIONS ********************/

/**
 * @brief Compute keypoints corresponding to the image then tries to match them to detect copy-move forgeries.
 *
 * @param w,h,channels: size of the image (respectively width, height and channels)
 * @param image: input image
 * @param data: used to return a condensed version of every matches
 * @param matchings: used to return the list of every matches 
 * @param ps: patch size used for the descriptors 
 * @param tau: the threshold for the matching, if automatic is true it corresponds to the number of false alarm
 * @param automatic: if true, the threshold is estmated using the NFA formula using tau as number of false alarm
 **/
double perform_matching(int channels, std::vector<float>& image, int w, int h, std::vector<float>& data, Matchingslist& matchings, int ps, float tau, bool automatic, bool verbose=false);


/**
 * @brief Compute matches between precomputed keypoints.
 *
 * @param c: number of channels of the image
 * @param keys: list of every precomputed keypoints
 * @param matchings: used to return the list of every matches 
 * @param ps: patch size used for the descriptors 
 * @param tau: the threshold for the matching, if automatic is true it corresponds to the number of false alarm
 * @param automatic: if true, the threshold is estmated using the NFA formula using tau as number of false alarm
 *
 * @return The total number of matches
 **/
double compute_matches(int c, std::vector<KeyPoint *> &keys, Matchingslist &matchings, int ps, float tau, bool automatic, bool verbose=false);

/**
 * @brief Computes all keypoints of a given image.
 *
 * @param image: input image.
 * @param width, height, channels: size of the image (respectively width, height and channels)
 * @param keypoints: list of keypoints that will be returned.
 * @param ps: patch size used for the descriptors 
 *
 * @return The total number of keypoints that have been found.
 */
int compute_keypoints(std::vector<float>& image, int width, int height, int channels, std::vector<KeyPoints*>& keypoints, int ps);

/**
 * @brief Compute a pair of descriptors.
 *
 * @param grad_*: gradients in both direction of the two descriptors
 * @param ps: patch size used for the descriptors 
 * @param ch: number of channels of the image
 * @param tau: the threshold for the matching
 * @param flip: if true test the flipped version of the descriptors
 *
 * @return The total number of matches
 **/
bool patch_comparison(double* grad_x_1, double* grad_y_1, double* grad_x_2, double* grad_y_2,
                         int ps, int ch, double tau, bool flip = false);


#endif // MATCHING_H
