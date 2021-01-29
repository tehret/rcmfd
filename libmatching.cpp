#include "matching.h"
#include "utils/drawing.h"

/**
 * @brief Creates an image showing the matches provided
 *
 * @param w,h,channels: size of the image (respectively width, height and channels)
 * @param im: image corresponding to the matches
 * @param matchings: list of matches to draw on the image
 **/
void show_images_matches(int channels, std::vector<float>& im, int ps, int w, int h, Matchingslist& matchings, float* rgb)
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
    for(int c=0;c<3;c++)
    for(int j = 0; j < (int) h; j++)
    for(int i = 0; i < (int) w; i++)
        rgb[j*w*3+i*3+c] = outimg[c][j*w+i];

    for(int c=0;c<3;c++)
        delete[] outimg[c];
}

#ifdef __cplusplus
extern "C"
#endif
bool wrapper_perform_matching(int c, float* image, int w, int h, int ps, float tau, bool automatic, float* output, bool verbose)
{
    Matchingslist matchings;
    std::vector<float> data;
    std::vector<float> im(image, image + h*w*c);

    double threshold = perform_matching(c, im, w, h, data, matchings, ps, tau, automatic, verbose);

    if(output != NULL)
        show_images_matches(c, im, ps, w, h, matchings, output);
    
    return (data.size() > 0);
}
