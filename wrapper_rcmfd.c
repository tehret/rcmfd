#include "libmatching.h"

bool perform_matching_py(int channels, float* image, int w, int h, int ps, float tau, bool automatic, float* output, bool verbose)
{
    return wrapper_perform_matching(channels, image, w, h, ps, tau, automatic, output, verbose);
}
