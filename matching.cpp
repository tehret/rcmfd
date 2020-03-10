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

#include "matching.h"

#define ABS(x)    (((x) > 0) ? (x) : (-(x)))

using namespace std;

void compute_local_descriptor_keypoints(vector<float>& queryImg, int w, int h, int c, Keypointlist& KPs, int ps)
{
    keypointslist* keys = new keypointslist;

    // Use fault SIFT parameters
    siftPar siftparameters;
    default_sift_parameters(siftparameters);
    NewOriSize1 = ps;

    // Compute keypoints
    compute_sift_keypoints(queryImg.data(), *keys, w, h, c, siftparameters);

    float* outimg = new float[w*h];
    for(int j = 0; j < (int) h; j++)
    for(int i = 0; i < (int) w; i++)
        outimg[j*w+i] = (c > 1) ? (queryImg[j*w+i] + queryImg[1*w*h+j*w+i] + queryImg[2*w*h+j*w+i])/3. : queryImg[j*w+i];

    // Save keypoints into a different structure
    KPs.resize(keys->size());
    for(int i=0; i<(int)keys->size();i++)
    {
        KPs[i].x = (*keys)[i].x;
        KPs[i].y = (*keys)[i].y;
        KPs[i].kp_ptr = &((*keys)[i]);
        KPs[i].scale = (*keys)[i].scale;
        KPs[i].angle = (*keys)[i].angle;

        int sq = (int)(KPs[i].scale * ps / 2);
        draw_square(outimg,  round(KPs[i].x)-sq, round(KPs[i].y)-sq, 2*sq, 2*sq, 255, w, h);
    }
}

void add_keypoint(Keypointlist& keys, std::vector<KeyPoints*>& mapKP, int width, int height)
{
    float x,y;
    int xr,yr;
    bool only_center;
    for(int i=0; i<(int)keys.size();i++)
    {
        x = keys[i].x;
        xr = (int) round(x);
        y = keys[i].y;
        yr = (int) round(y);

        int ind =  yr*width + xr, newind;
        newind = ind;
        // Order matches by their position for efficient access
        if ( mapKP[ind]==0 )
        {
            mapKP[ind] = new KeyPoints();
            mapKP[ind]->x = x;
            mapKP[ind]->y = y;
            mapKP[ind]->sum_x = x;
            mapKP[ind]->sum_y = y;
            mapKP[ind]->KPvec.push_back(keys[i]);
        }
        else
        {
            mapKP[ind]->KPvec.push_back(keys[i]);
            mapKP[ind]->sum_x += x;
            mapKP[ind]->sum_y += y;
            mapKP[ind]->x = mapKP[ind]->sum_x / mapKP[ind]->KPvec.size();
            mapKP[ind]->y = mapKP[ind]->sum_y / mapKP[ind]->KPvec.size();

        }

        only_center = false;
        int r = FILTER_RADIUS;
        while(!only_center)
        {
            only_center = true;
            for (int xi = xr-r; xi<=xr+r; xi++)
                for (int yi = yr-r; yi<=yr+r; yi++)
                {
                    if (( sqrt(pow(xi-xr,2) + pow(yi-yr,2))<=r )&&(xi>0)&&(xi<width)&&(yi>0)&&(yi<height))
                    {
                        int indi = yi*width + xi;
                        // Group keypoints that are very close to each other
                        if ( (mapKP[indi]!=0)&&(indi!=ind) )
                        {
                            // Merge indi to ind
                            only_center = false;
                            mapKP[ind]->sum_x += mapKP[indi]->sum_x;
                            mapKP[ind]->sum_y += mapKP[indi]->sum_y;

                            for (int k=0;k<(int)mapKP[indi]->KPvec.size();k++)
                            {
                                mapKP[ind]->KPvec.push_back(mapKP[indi]->KPvec[k]);
                            }
                            mapKP[indi] = 0;

                            mapKP[ind]->x = mapKP[ind]->sum_x / mapKP[ind]->KPvec.size();
                            mapKP[ind]->y = mapKP[ind]->sum_y / mapKP[ind]->KPvec.size();

                            // Update newind
                            x = mapKP[ind]->x;
                            xr = (int) round(x);
                            y = mapKP[ind]->y;
                            yr = (int) round(y);
                            newind = yr*width + xr;

                            if (newind!=ind)
                            {
                                if (mapKP[newind]==0)
                                { // newind empty
                                    mapKP[newind] = mapKP[ind];
                                    mapKP[ind] = 0;
                                    ind = newind;
                                }
                                else
                                { // newind not empty
                                    mapKP[newind]->sum_x += mapKP[ind]->sum_x;
                                    mapKP[newind]->sum_y += mapKP[ind]->sum_y;

                                    for (int k=0;k<(int)mapKP[ind]->KPvec.size();k++)
                                    {
                                        mapKP[newind]->KPvec.push_back(mapKP[ind]->KPvec[k]);
                                    }

                                    mapKP[newind]->x = mapKP[newind]->sum_x / mapKP[newind]->KPvec.size();
                                    mapKP[newind]->y = mapKP[newind]->sum_y / mapKP[newind]->KPvec.size();
                                    mapKP[ind] = 0;
                                    ind = newind;
                                }
                            }
                        }
                    }

                }
        }

    }

}


int compute_keypoints(vector<float>& image, int width, int height, int channels, std::vector<KeyPoints*>& keypoints, int ps)
{	

    std::vector<KeyPoints*> mapKP;
    mapKP.resize(width*height);

    for(int i=0;i<width*height;i++)
        mapKP[0]= 0;

    int num_keys_total=0;

    Keypointlist keys;
    compute_local_descriptor_keypoints(image, width, height, channels, keys, ps);
    add_keypoint(keys, mapKP, width, height);

    // save in `keypoints`
    int num_max = 0, num_min = 500000, total = 0;
    float num_mean = 0;
    for (int i = 0; i < (int) mapKP.size(); i++)
        if (mapKP[i]!=0)
        {
            keypoints.push_back(mapKP[i]);
            total +=mapKP[i]->KPvec.size();
            num_keys_total += 1;
            if (num_max<(int)mapKP[i]->KPvec.size())
                num_max = mapKP[i]->KPvec.size();
            if (num_min>(int)mapKP[i]->KPvec.size())
                num_min = mapKP[i]->KPvec.size();
            num_mean +=mapKP[i]->KPvec.size();
        }
    num_mean = num_mean/num_keys_total;

    return num_keys_total;
}

bool patch_comparison(double * grad_x_1, double * grad_y_1, double * grad_x_2, double * grad_y_2,
        int ps, int ch, double tau, bool flip )
{
    for(int c = 0; c < ch; ++c)
        for(int y = 1; y < ps-1; ++y)
        {
            for(int x = 1; x < ps-1; ++x)
            {
                double dist;
                // Check whether patches are well defined
                if(grad_x_1[x+y*ps+c*ps*ps] == NOTDEF || grad_x_2[x+y*ps+c*ps*ps] == NOTDEF || grad_y_1[x+y*ps+c*ps*ps] == NOTDEF || grad_y_2[x+y*ps+c*ps*ps] == NOTDEF)
                    return false;
                if(flip)
                    dist = (grad_x_1[x+(ps-y-1)*ps+c*ps*ps] - grad_x_2[x+y*ps+c*ps*ps])*(grad_x_1[x+(ps-y-1)*ps+c*ps*ps] - grad_x_2[x+y*ps+c*ps*ps]) + 
                        (grad_y_1[x+(ps-y-1)*ps+c*ps*ps] + grad_y_2[x+y*ps+c*ps*ps])*(grad_y_1[x+(ps-y-1)*ps+c*ps*ps] + grad_y_2[x+y*ps+c*ps*ps]);
                else
                    dist = (grad_x_1[x+y*ps+c*ps*ps] - grad_x_2[x+y*ps+c*ps*ps])*(grad_x_1[x+y*ps+c*ps*ps] - grad_x_2[x+y*ps+c*ps*ps]) + 
                        (grad_y_1[x+y*ps+c*ps*ps] - grad_y_2[x+y*ps+c*ps*ps])*(grad_y_1[x+y*ps+c*ps*ps] - grad_y_2[x+y*ps+c*ps*ps]);

                if(dist > tau)
                    return false;
            }
        }
    return true;
}

double compute_matches(int c, std::vector<KeyPoints*>& keys, Matchingslist &matchings, int ps, float epsilon, bool automatic)
{  	
    int tstart = time(0);
    printf("Keypoints matching...\n");

    double tau;
    // Automatic estimation of the threshold if necessary
    if(automatic)
        tau = 2*stats::qchisq(pow(exp(epsilon)/(keys.size()*keys.size()), 1./((ps-2)*(ps-2)*c)), 1);
    else
        tau = epsilon;

    // Test every pair of descriptors
    for (int idx1 = 0; idx1 < (int)keys.size(); ++idx1)
    for (int idx2 = 0; idx2 < idx1; ++idx2)
    {
        // Avoid matching a descriptor with itself
        if(idx1 == idx2)
            continue;

        bool stopped = false;
        for(int i1 = 0; i1 < (int)keys[idx1]->KPvec.size() && !stopped; ++i1)
        for(int i2 = 0; i2 < (int)keys[idx2]->KPvec.size(); ++i2)
        {
            // Check overlapping descriptors
            if(std::sqrt((keys[idx1]->x - keys[idx2]->x)*(keys[idx1]->x - keys[idx2]->x) + (keys[idx1]->y - keys[idx2]->y)*(keys[idx1]->y - keys[idx2]->y)) < ps/2*(keys[idx1]->KPvec[i1].scale + keys[idx2]->KPvec[i2].scale))
                continue;

            double sigma2 = 1;
            if(automatic)
            {
               sigma2 = std::min(sigma2,std::min(static_cast<keypoint*>(keys[idx1]->KPvec[i1].kp_ptr)->var, static_cast<keypoint*>(keys[idx2]->KPvec[i2].kp_ptr)->var)); 
            }

            // Try to match both descriptors and when one if flipped 
            if(patch_comparison(
                        static_cast<keypoint*>(keys[idx1]->KPvec[i1].kp_ptr)->gradx,
                        static_cast<keypoint*>(keys[idx1]->KPvec[i1].kp_ptr)->grady,
                        static_cast<keypoint*>(keys[idx2]->KPvec[i2].kp_ptr)->gradx,
                        static_cast<keypoint*>(keys[idx2]->KPvec[i2].kp_ptr)->grady,
                        ps, c, sigma2*tau, false) ||
                    patch_comparison(
                        static_cast<keypoint*>(keys[idx1]->KPvec[i1].kp_ptr)->gradx,
                        static_cast<keypoint*>(keys[idx1]->KPvec[i1].kp_ptr)->grady,
                        static_cast<keypoint*>(keys[idx2]->KPvec[i2].kp_ptr)->gradx,
                        static_cast<keypoint*>(keys[idx2]->KPvec[i2].kp_ptr)->grady,
                        ps, c, sigma2*tau, true))
            {
                Keypoint_simple k1, k2;

                k1.x = keys[idx1]->x;
                k1.y = keys[idx1]->y;
                k1.scale = keys[idx1]->KPvec[i1].scale;
                k1.angle = keys[idx1]->KPvec[i1].angle;

                k2.x = keys[idx2]->x;
                k2.y = keys[idx2]->y;
                k2.scale = keys[idx2]->KPvec[i2].scale;
                k2.angle = keys[idx2]->KPvec[i2].angle;

                matchings.push_back( Matching(k1,k2) );
                stopped = true;
                break;
            }
        }
    }

    printf("   %d matches found.\n", (int) matchings.size());
    printf("Matching accomplished in %ld seconds.\n\n", (time(0) - tstart));
    return tau;
}


double perform_matching(int channels, vector<float>& image, int w, int h, vector<float>& data, Matchingslist& matchings, int ps, float tau, bool automatic)
{
    // Compute keypoints
    std::vector<KeyPoints*> keys;
    keys.clear();

    int num_keys=0;

    int tstart = time(0);
    num_keys = compute_keypoints(image, w, h, channels, keys, ps);
    printf("%d keypoints have been found\n", num_keys);
    printf("Keypoints computation accomplished in %ld seconds.\n \n", time(0) - tstart);


    // Match keypoints
    double threshold = compute_matches(channels, keys, matchings, ps, tau, automatic);

    // Generate data matrix: there are eight rows of info
    for ( int i = 0; i < (int) matchings.size(); i++ )
    {
        Matching *ptr_in = &(matchings[i]);
        data.push_back(ptr_in->first.x);
        data.push_back(ptr_in->first.y);
        data.push_back(ptr_in->first.scale);
        data.push_back(ptr_in->first.angle);
        data.push_back(ptr_in->second.x);
        data.push_back(ptr_in->second.y);
        data.push_back(ptr_in->second.scale);
        data.push_back(ptr_in->second.angle);
    }
    printf("Done.\n\n");

    if(matchings.size() > 0)
        printf("This image IS forged\n");
    else
        printf("This image IS NOT forged\n");

    keys.clear();
    return threshold;
}
