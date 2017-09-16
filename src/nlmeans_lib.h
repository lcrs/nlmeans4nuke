/* nlmeans++_lib.cpp */

#ifndef __cplusplus
#error "a C++ compiler is expected"
#endif

using namespace std;
#include <vector>
#include <list>
#include <DDImage/Iop.h>

void noise_estimate_compute_derivative(float *igray, float *igray2, int width, int height);
void noise_estimate_derivative_variance(float *igray, float *ider, int s,   list<float> &stdvalues, list<float> &averages, int width, int height);
void noise_estimation_algorithm(int iBins, int iWsize, float percent, int equiflag,float *xresults, float *results, float *igray, int width, int height);
void nlmeans(int nwin,int bloc,int averageoption, float *kernel, int kwidth, float multiplier, float rsigma,float *ired,float * ored,int width,int height);
void nlmeans(int nwin, int bloc, int averageoption, float *kernel, int ksize, float multiplier, float rsigma, float gsigma, float bsigma,float *ired,float *igreen,float * iblue,float * ored,float * ogreen,float * oblue,int width,int height, DD::Image::Iop *myop);
void nlmeans_sliding_window(int nwin,int bloc, float *kernel, int kwidth, float multiplier, float rsigma, float gsigma, float bsigma,float *ired,float *igreen,float * iblue,float * ored,float * ogreen,float * oblue,int width,int height);

