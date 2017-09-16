/*
 * Copyright 2009, 2010 IPOL Image Processing On Line http://www.ipol.im/
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file nlmeans_lib.cpp
 * @brief denoising functions
 *
 * @author Antoni Buades <toni.buades@uib.es>
 */

/*

Modifed slightly 08/2010 for use in Nuke plugin by lewis@lewissaunders.com.

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "nlmeans_lib.h"

#define MAX(i,j) ( (i)<(j) ? (j):(i) )
#define MIN(i,j) ( (i)<(j) ? (i):(j) )

#define LUTMAX 30
#define LUTPRECISION 1000.0

#define TINY 1.0e-10
#define MAXFLOAT 10000000.0

static void clear(float *u,float value,int size) { for(int i=0; i<size; i++) u[i]=value;  }

/* Quicksort,  values in arr are set in increasing order and brr elements are switched at the same time*/
static void FSWAP(float *x,float *y)
{
  float aux;
  aux=*x;
  *x=*y;
  *y=aux;
}



static void quick_sort(float *arr,float *brr,int n)
{
  int M=7,NSTACK=50;
  int i,ir,j,k,jstack=-1,l=0;
  float a,b;
  int istack[50];
  
  ir=n-1;


  for(;;){
    if(ir-l<M){
      for(j=l+1;j<=ir;j++){
	a=arr[j];
	b=brr[j];
	for(i=j-1;i>=l;i--){
	  if (arr[i]<=a) break;
	  arr[i+1]=arr[i];
	  brr[i+1]=brr[i];
	}
	arr[i+1]=a;
	brr[i+1]=b;

      }

      if (jstack<0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {

      k=(l+ir) >> 1;
      FSWAP(&arr[k],&arr[l+1]);
      FSWAP(&brr[k],&brr[l+1]);
      if (arr[l]>arr[ir]){
	FSWAP(&arr[l],&arr[ir]);
	FSWAP(&brr[l],&brr[ir]);
      }
      if (arr[l+1]>arr[ir]){
	FSWAP(&arr[l+1],&arr[ir]);
	FSWAP(&brr[l+1],&brr[ir]);
      }
      if (arr[l]>arr[l+1]){
	FSWAP(&arr[l],&arr[l+1]);
	FSWAP(&brr[l],&brr[l+1]);
      }
      i=l+1;
      j=ir;
      a=arr[l+1];
      b=brr[l+1];
      for(;;){
	do i++; while (arr[i]<a);
	do j--; while (arr[j]>a);
	if (j<i) break;
	FSWAP(&arr[i],&arr[j]);
	FSWAP(&brr[i],&brr[j]);
      }

      arr[l+1]=arr[j];
      arr[j]=a;
      brr[l+1]=brr[j];
      brr[j]=b;
      jstack+=2;

      if (jstack>=NSTACK) { printf("Stack too small\n"); exit(-1);}
      if (ir-i+1>=j-l){
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }

    }
  }

  
}




void noise_estimate_derivative_variance(float *igray, float *ider, int s,  list<float> & stdvalues, list<float> & averages, int width, int height)
{

	int blockSize = s*s;
	
	for(int y=0; y + s < height; y++)
		for(int x=0; x + s < width; x++)
	{

		////////////////// Computing std for each block
		float dmean = 0.0;
		float dmean2 = 0.0;
		float mean = 0.0;

		for(int i = 0; i < s; i++)
			for(int j = 0; j < s; j++)
		{
	
			float value = igray[(y+j)*width + (x+i)];
			mean += value;

			float dvalue = ider[(y+j)*width + (x+i)];

			dmean += dvalue;
			dmean2 += dvalue * dvalue;

		}

		mean /= (float) blockSize;
		dmean /= (float) blockSize;
		dmean2 /= (float) blockSize;

		float std = sqrtf(dmean2 - dmean*dmean);

		//////////////////
		stdvalues.push_back(std);	
		averages.push_back(mean);			
	}
}




void noise_estimate_compute_derivative(float *igray, float *igray2, int width, int height)
{

	float sqrt2 = sqrtf(2.0);
	for(int y=1; y < height - 1; y++)
		for(int x=1; x < width - 1; x++)
	{

		////////////////// Computing derivative for each point
// 		float derivative = 4.0*igray[y*width+x] - igray[y*width+x-1] - igray[(y+1)*width+x] - igray[y*width+x+1] - igray[(y-1)*width+x] ;	
// 		derivative /= sqrt20;

		float derivative = igray[y*width+x] - igray[y*width+x+1];	
		derivative /= sqrt2;

		igray2[y*width+x] = derivative;

	}
}







void noise_estimation_algorithm(int iBins, int iWsize, float percent, int equiflag,float *xresults, float *results, float *igray, int width, int height)
{

	float correction = 1.0;
	if (iWsize == 15 && percent == 0.01f) correction = 1.1502;
	else if (iWsize == 21 && percent == 0.01f) correction = 1.1035;
	else if (iWsize == 15 && percent == 0.005f) correction = 1.1676;
	else if (iWsize == 21 && percent == 0.005f) correction = 1.1153;
	else {  printf("ERROR :: noise_estimation_algorithm :: parameters not calibrated\n"); exit(-1); }

	//// Computing derivatives
	float *ider = new float[width*height];
	clear(ider,0.0,width*height); 
	noise_estimate_compute_derivative(igray,ider,width,height);

	
	///// Partial standard deviation values and averages
	list<float> stdvalues;
	list<float> averages;
	stdvalues.clear(); averages.clear();

	///// Compute variance of patches of derivatives
	noise_estimate_derivative_variance(igray, ider,  iWsize,   stdvalues, averages, width, height);



	////  Order averages::  averages --> orderaverages
	int asize = averages.size();
	float* order_averages = new float[asize];
	float* order_stdvalues = new float[asize];

	list<float>::iterator aptr = averages.begin();
	list<float>::iterator sptr = stdvalues.begin();
 	for(int j=0; j < asize; aptr++, sptr++,j++) 
 	{
 		order_averages[j] = *aptr; 
 		order_stdvalues[j] =  *sptr;
 	}


	quick_sort(order_averages, order_stdvalues, asize);


	
	//// Limits of quantization steps
	float * hlimits = new float[iBins+1];	
	clear(hlimits, 0.0, iBins+1);
	
	float rMinValue = order_averages[0];
	float rMaxValue = order_averages[asize-1]; 
		 
		
	int cStep;
	if (equiflag)
	{
		int uasize = 0;
		for(; order_averages[uasize] < rMaxValue; uasize++) {;}

		//// Step histogram quantization
		int iStep = (int) ceilf((float) (uasize) / (float) (iBins));
		cStep = iStep;

		int k=0;
		for(int i = 0; i < uasize; i++)
			if (i % iStep == 0 && k < iBins)
			{
				hlimits[k] = order_averages[i];
				k++;
			}

		hlimits[iBins] = rMaxValue;
	}	
	else
	{


		float rStep = (rMaxValue - rMinValue) / (float) iBins;	
		cStep = (int) rStep;

		for(int i=0; i < iBins; i++)	hlimits[i] = rMinValue + (float) i * rStep;
		hlimits[iBins] = rMaxValue;

	}


	//// order std inside each range bin, means ordered according to std
	float * partial_order_averages = new float[asize];
	float * partial_order_stdvalues = new float[asize];

	for(int i=0, k=0; i < iBins; i++)
	{

		int counter = 0;
		while(order_averages[k] >= hlimits[i] && order_averages[k] < hlimits[i+1])
		{
			partial_order_averages[counter] = order_averages[k];
			partial_order_stdvalues[counter] = order_stdvalues[k];
			counter++;

			k++;
		}

		quick_sort(partial_order_stdvalues, partial_order_averages, counter);


		int npercentile = (int) (percent * (float) counter);
		
		results[i] = correction * partial_order_stdvalues[npercentile];
		xresults[i] =  partial_order_averages[npercentile];
			
	}	


// 	for(int i=0, k=0; i < iBins; i++) printf("%d: %f \t",i,results[i]);
// 	printf("\n");

}

static int normalize(float *u,int size) 
{
  
  float some = 0.0;  
  for(int i=0;i<size;i++) some+=u[i];
  if (some != 0.0) {	for(int i=0;i<size;i++) u[i]/=some;  return 1;}  
  else return 0;	
}



static float l2_distance(float *u0,float *u1,int i0,int j0,int i1,int j1,int radius,int width)
{

	int wsize=(2*radius+1)*(2*radius+1);

	float dist=0.0;       
	for (int s=-radius; s<= radius; s++){
	
		int l = (j0+s)*width + (i0-radius);
		float *ptr0 = &u0[l];
	
		l = (j1+s)*width + (i1-radius);
		float *ptr1 = &u1[l];
	
		for(int r=-radius;r<=radius;r++,ptr0++,ptr1++){	float dif = (*ptr0 - *ptr1); dist += (dif*dif); }

	}

	dist/=(float) wsize;
	return dist;
}


static float weighted_l2_distance(float *u0,float *u1,int i0,int j0,int i1,int j1,int width,float *kernel,int radius)
{


	float *ptrk=&kernel[0];
	float dist=0.0;       
	for (int s=-radius; s<= radius; s++){
	
		int l = (j0+s)*width + (i0-radius);
		float *ptr0 = &u0[l];
	
		l = (j1+s)*width + (i1-radius);
		float *ptr1 = &u1[l];
	
	
		for(int r=-radius;r<=radius;r++,ptr0++,ptr1++,ptrk++){ float dif = (*ptr0 - *ptr1); dist += *ptrk*(dif*dif); }
	
	}

	return dist;
}






/////   research zone: 2*bloc+1 x 2*bloc+1
////    comparison window: 	if kernel == NULL 	flat of size 2*nwin+1 x 2*nwin+1 
////				else	
////					weighted by kernel and size kwidth*kwidth
////
////	rsigma, gsigma, bsigma: noise standard deviation for each channel 

void nlmeans(int nwin,int bloc,int averageoption, float *kernel, int kwidth, float multiplier, float rsigma,float *ired,float * ored,int width,int height)
{
	
	int wxh = width*height;


	/////// Clear output
	clear(ored,0.0,wxh);

	/////// Normalize kernel if necessary
	if (kernel){
	   
		nwin = (kwidth-1) / 2;    
		normalize(kernel,kwidth*kwidth);
	}




	///// Mean of squared distances of patches centered in (x0,y0) with added noise of variance sigma equals 2 sigma ^2 + d(x0,y0)
	///// Variance of squared distance of patches centered in (x0,y0) with added noise of variance sigma equals (8 sigma^4 / N ) + ( 8 sigma^2 d (x0,y0) / N)
	/////		- N number of pixels 27 (3x3x3), 75 (3x5x5)
	///// 		- d(x0,y0) distance of non-noisy patches 
 	float sigma2 = rsigma*rsigma;
	float fnpix = (float) ((2*nwin+1) * (2*nwin+1));
 	float dx0y0 = 0.0f;


	float variance=0.0;
	if (kernel)
	{
		for(int i=0; i < kwidth*kwidth; i++) variance += kernel[i]*kernel[i];
		variance *= 8.0f * sigma2* sigma2;

	} else
 		variance = (8.0f * sigma2 * sigma2 / fnpix) + (8.0f * sigma2 * dx0y0 / fnpix);
 
	

	float std = sqrtf(variance);
 	float t1 = 2.0f * sigma2 ;  				// mean
	float t2 = 2.0f * sigma2 + multiplier * std;  		// mean + 2.0 std
	float t2mt1 = t2 - t1 ;					// useful for weight computation

	
	////////// We denoise the entire comparison window and not only the central pixel
	float red = 0.;
	float *rwindow = NULL;
	int winsize = (2*nwin+1)*(2*nwin+1);
	if (averageoption)
	{
		rwindow = new float[winsize];
	}


	/////////  In general each estimate will have nine possible values except for boundary pixels
	////////   We carry for each pixel the number of estimates we dispose 
	float *counter_mask = new float[wxh];
	clear(counter_mask,0.0,wxh);


	//////////////  PROCESS STARTS: for each pixel (x,y)
	for(int y=0; y < height ; y++)
		for(int x=0 ; x < width; x++){
      

		
		////////// We reduce the size of the comparison window if we are near the boundary  
		int nwin0 = MIN(nwin,MIN(width-1-x,MIN(height-1-y,MIN(x,y))));

  
		////////// Learning zone depending on the boundary and the size of the window
		int imin=MAX(x-bloc,nwin0);
		int jmin=MAX(y-bloc,nwin0);

		int imax=MIN(x+bloc,width-1-nwin0);
		int jmax=MIN(y+bloc,height-1-nwin0);
 
      
		/////////// Clear current denoised window
		if (averageoption)
		{
			clear(rwindow,0.0,winsize);
		} else
		{
			red = 0.0; 
		}


		float totalweight=0.0;

		for(int j=jmin; j <= jmax; j++)
			for(int i=imin ; i <= imax; i++) {
		
				float dif = 0.0f;


				if (kernel && nwin0==nwin){

					dif = weighted_l2_distance(ired,ired,x,y,i,j,width,kernel,nwin);
		
				} else
				{
					dif=l2_distance(ired,ired,x,y,i,j,nwin0,width);
				}

		
				/////// If difference less than threshold average window (i,j)
				if (dif < t2)	
				{					

					float weight;

					///// weight => 1 		if dist < t1
					////		(t2-x)/(t2-t1)	if t1<dist<t2
					////		0		if t>t2
				
					//// weight function is continuous but not derivable at t1 and t2
					 
					if (dif < t1) weight = 1.0f;
					else weight = (t2 - dif) / t2mt1;	// t2mt1 = t2 - t1 (declared outside fors)


					if (averageoption) { 

						///// We fill window of size 2*nwin0+1 x 2*nwin0+1
						for(int r=-nwin0; r<=nwin0; r++)
							for(int s=-nwin0; s <= nwin0; s++)
						{ 
	
							int index = (nwin+s) * (2*nwin+1) + r + nwin;	//// However the complete rwindow/gwindow/bwindow have size depending on nwin
							int l=(j+s)*width+i+r;
	
							rwindow[index] +=   weight * ired[l];
		
						}		

					} else {

							int l=j*width+i;	
							red +=   weight * ired[l];
					}


					//// OJO ///// Compute totalweight for each pixel and normalize finally for each pixel --> Should be similar to Stationary NLmeans
					totalweight += weight;		
				}


			}
		  
            
			if (totalweight > TINY) {
			
				if (averageoption) { 

				     for(int r=-nwin0; r <= nwin0; r++)
					for(int s=-nwin0; s <= nwin0;s++)
					{ 
						int index = (nwin+s) *  (2*nwin+1) +r+nwin;
						int l=(y+s)*width+x+r;

						//// Divide estimate by the sum of weights and carry this value for this pixel
						ored[l] += (rwindow[index] / totalweight); 

						counter_mask[l]++;
					}

				} else {

					int l=y*width+x;

					//// Divide estimate by the sum of weights and carry this value for this pixel
					ored[l] = (red / totalweight); 
 
				}	

			}

      
		}


		if (averageoption){


			for(int i=0; i < wxh;i++)
				if (counter_mask[i]>0.0) 
				{
					ored[i] /= counter_mask[i];
	
				} else  //// This should never happen, since at least for each pixel itself has distance 0
				{
					ored[i] = ired[i];
				}

		}


		////////// We put back a 15% of the original image
		for(int i=0; i < wxh;i++)
		{
			ored[i] = 0.85f * ored[i] + 0.15f * ired[i];
		}


}






/////   research zone: 2*bloc+1 x 2*bloc+1
////    comparison window: 	if kernel == NULL 	flat of size 2*nwin+1 x 2*nwin+1 
////				else	
////					weighted by kernel and size kwidth*kwidth
////
////	rsigma, gsigma, bsigma: noise standard deviation for each channel 

void nlmeans(int nwin,int bloc,int averageoption, float *kernel, int kwidth, float multiplier, float rsigma, float gsigma, float bsigma,float *ired,float *igreen,float * iblue,float * ored,float * ogreen,float * oblue,int width,int height, DD::Image::Iop *myop)
{
	
	int wxh = width*height;


	/////// Clear output
	clear(ored,0.0,wxh);
	clear(ogreen,0.0,wxh);
	clear(oblue,0.0,wxh);

	/////// Normalize kernel if necessary
	if (kernel){
	   
		nwin = (kwidth-1) / 2;    
		normalize(kernel,kwidth*kwidth);
	}


	////// Noise details and information
	float sigma2 = (rsigma*rsigma + bsigma*bsigma + gsigma*gsigma) / 3.0f;


	///// Mean of squared distances of patches centered in (x0,y0) with added noise of variance sigma equals 2 sigma ^2 + d(x0,y0)
	///// Variance of squared distance of patches centered in (x0,y0) with added noise of variance sigma equals (8 sigma^4 / N ) + ( 8 sigma^2 d (x0,y0) / N)
	/////		- N number of pixels 27 (3x3x3), 75 (3x5x5)
	///// 		- d(x0,y0) distance of non-noisy patches 
 	float fnpix = (float) (3 * (2*nwin+1) * (2*nwin+1));
 	float dx0y0 = 0.0f;


	float variance=0.0;
	if (kernel)
	{
		for(int i=0; i < kwidth*kwidth; i++) variance += kernel[i]*kernel[i];
		
		variance = variance * 8.0f * sigma2 * sigma2 / 9.0f;
	
	} else
 		variance = (8.0f * sigma2 * sigma2 / fnpix) + (8.0f * sigma2 * dx0y0 / fnpix);
 
 	float std = sqrtf(variance);
 	float t1 = 2.0f * sigma2 ;  				// mean
	float t2 = 2.0f * sigma2 + multiplier * std;  		// mean + 2.0 std
	float t2mt1 = t2 - t1 ;					// useful for weight computation
	
	////////// We denoise the entire comparison window and not only the central pixel
	float red = 0., green = 0., blue = 0.;
	float *rwindow = NULL, *gwindow = NULL, *bwindow = NULL;
	int winsize = (2*nwin+1)*(2*nwin+1);
	if (averageoption)
	{
		rwindow = new float[winsize];
		gwindow = new float[winsize];
		bwindow = new float[winsize];
	}


	/////////  In general each estimate will have nine possible values except for boundary pixels
	////////   We carry for each pixel the number of estimates we dispose 
	float *counter_mask = new float[wxh];
	clear(counter_mask,0.0,wxh);


	//////////////  PROCESS STARTS: for each pixel (x,y)
	for(int y=0; y < height; y++) {
		myop->progressFraction((double)y / (double)height);
		if(myop->aborted()) return;
		
		for(int x=0 ; x < width; x++){
      

		
		////////// We reduce the size of the comparison window if we are near the boundary  
		int nwin0 = MIN(nwin,MIN(width-1-x,MIN(height-1-y,MIN(x,y))));

  
		////////// Learning zone depending on the boundary and the size of the window
		int imin=MAX(x-bloc,nwin0);
		int jmin=MAX(y-bloc,nwin0);

		int imax=MIN(x+bloc,width-1-nwin0);
		int jmax=MIN(y+bloc,height-1-nwin0);
 
      
		/////////// Clear current denoised window
		if (averageoption)
		{
			clear(rwindow,0.0,winsize);
			clear(gwindow,0.0,winsize);
			clear(bwindow,0.0,winsize);
		} else
		{
			red = 0.0; green = 0.0; blue = 0.0;
		}


		float totalweight=0.0;

		for(int j=jmin; j <= jmax; j++)
			for(int i=imin ; i <= imax; i++) {
		
				float dif = 0.0f;


				if (kernel && nwin0==nwin){

					dif = weighted_l2_distance(ired,ired,x,y,i,j,width,kernel,nwin);
					dif += weighted_l2_distance(igreen,igreen,x,y,i,j,width,kernel,nwin);
					dif += weighted_l2_distance(iblue,iblue,x,y,i,j,width,kernel,nwin);
		

				} else
				{

					dif  =  l2_distance(igreen,igreen,x,y,i,j,nwin0,width);
					dif += l2_distance(ired,ired,x,y,i,j,nwin0,width);
					dif += l2_distance(iblue,iblue,x,y,i,j,nwin0,width);
					
				}

				dif = dif / 3.0;
				/////// If difference less than threshold average window (i,j)
				if (dif < t2)	
				{					

					float weight;

					///// weight => 1 		if dist < t1
					////		(t2-x)/(t2-t1)	if t1<dist<t2
					////		0		if t>t2
				
					//// weight function is continuous but not derivable at t1 and t2
					 
					if (dif < t1) weight = 1.0f;
					else weight = (t2 - dif) / t2mt1;	// t2mt1 = t2 - t1 (declared outside fors)


					if (averageoption) { 

						///// We fill window of size 2*nwin0+1 x 2*nwin0+1
						for(int r=-nwin0; r<=nwin0; r++)
							for(int s=-nwin0; s <= nwin0; s++)
						{ 
	
							int index = (nwin+s) * (2*nwin+1) + r + nwin;	//// However the complete rwindow/gwindow/bwindow have size depending on nwin
							int l=(j+s)*width+i+r;
	
							rwindow[index] +=   weight * ired[l];
							gwindow[index] +=   weight * igreen[l];
							bwindow[index] +=   weight * iblue[l];	
	
						}		

					} else {

							int l=j*width+i;
	
							red +=   weight * ired[l];
							green +=   weight * igreen[l];
							blue +=   weight * iblue[l];	
					}


					//// OJO ///// Compute totalweight for each pixel and normalize finally for each pixel --> Should be similar to Stationary NLmeans
					totalweight += weight;		
				}


			}
		  
            
			if (totalweight > TINY) {
			
				if (averageoption) { 

				     for(int r=-nwin0; r <= nwin0; r++)
					for(int s=-nwin0; s <= nwin0;s++)
					{ 
						int index = (nwin+s) *  (2*nwin+1) +r+nwin;
						int l=(y+s)*width+x+r;

						//// Divide estimate by the sum of weights and carry this value for this pixel
						ored[l] += (rwindow[index] / totalweight); 
						ogreen[l] += (gwindow[index] / totalweight); 
						oblue[l] += (bwindow[index] / totalweight); 

						counter_mask[l]++;
					}

				} else {

					int l=y*width+x;

					//// Divide estimate by the sum of weights and carry this value for this pixel
					ored[l] = (red / totalweight); 
					ogreen[l] = (green / totalweight); 
					oblue[l] = (blue / totalweight);
 
				}	

			}

      
		}

	}
		if (averageoption){


			for(int i=0; i < wxh;i++)
				if (counter_mask[i]>0.0) 
				{
					ored[i] /= counter_mask[i];
					ogreen[i] /= counter_mask[i];
					oblue[i] /= counter_mask[i];
	
				} else  //// This should never happen, since at least for each pixel itself has distance 0
				{
					ored[i] = ired[i];
					ogreen[i] = igreen[i];
					oblue[i] = iblue[i];
				}

		}
	

		/*////////// We put back a 15% of the original image
		for(int i=0; i < wxh;i++)
		{
			ored[i] = 0.85f * ored[i] + 0.15f * ired[i];
			ogreen[i] = 0.85f * ogreen[i] + 0.15f * igreen[i];
			oblue[i] = 0.85f * oblue[i] + 0.15f * iblue[i];	
		}*/

}





void nlmeans_sliding_window(int nwin,int bloc, float *kernel, int kwidth, float multiplier, float rsigma, float gsigma, float bsigma,float *ired,float *igreen,float * iblue,float * ored,float * ogreen,float * oblue,int width,int height)
{
	
	int wxh = width*height;


	/////// Clear output
	clear(ored,0.0,wxh);
	clear(ogreen,0.0,wxh);
	clear(oblue,0.0,wxh);

	/////// Normalize kernel if necessary
	if (kernel){
	   
		nwin = (kwidth-1) / 2;    
		normalize(kernel,kwidth*kwidth);
	}


	////// Noise details and information
	float sigma2 = (rsigma*rsigma + bsigma*bsigma + gsigma*gsigma) / 3.0f;


	///// Mean of squared distances of patches centered in (x0,y0) with added noise of variance sigma equals 2 sigma ^2 + d(x0,y0)
	///// Variance of squared distance of patches centered in (x0,y0) with added noise of variance sigma equals (8 sigma^4 / N ) + ( 8 sigma^2 d (x0,y0) / N)
	/////		- N number of pixels 27 (3x3x3), 75 (3x5x5)
	///// 		- d(x0,y0) distance of non-noisy patches 
 	float fnpix = (float) (3 * (2*nwin+1) * (2*nwin+1));
 	float dx0y0 = 0.0f;
 	float variance = (8.0f * sigma2 * sigma2 / fnpix) + (8.0f * sigma2 * dx0y0 / fnpix);
 	float std = sqrtf(variance);
 	float t1 = 2.0f * sigma2 ;  				// mean
	float t2 = 2.0f * sigma2 + multiplier * std;  		// mean + 2.0 std
	float t2mt1 = t2 - t1 ;					// useful for weight computation
	
// 	////////// We denoise the entire comparison window and not only the central pixel
// 	float red, green, blue;
// 	float *rwindow, *gwindow, *bwindow;
// 	int winsize = (2*nwin+1)*(2*nwin+1);
// 
// 	rwindow = new float[winsize];
// 	gwindow = new float[winsize];
// 	bwindow = new float[winsize];
//  


	/////////  In general each estimate will have nine possible values except for boundary pixels
	////////   We carry for each pixel the number of estimates we dispose 
	float *counter_mask = new float[wxh];
	clear(counter_mask,0.0,wxh);


	//////////////  PROCESS STARTS: for each pixel (x,y)
	for(int y=0; y < height ; y++)
		for(int x=0 ; x < width; x++){
      

		
		////////// We reduce the size of the comparison window if we are near the boundary  
		int nwin0 = MIN(nwin,MIN(width-1-x,MIN(height-1-y,MIN(x,y))));

  
		////////// Learning zone depending on the boundary and the size of the window
		int imin=MAX(x-bloc,nwin0);
		int jmin=MAX(y-bloc,nwin0);

		int imax=MIN(x+bloc,width-1-nwin0);
		int jmax=MIN(y+bloc,height-1-nwin0);
 
      /*
		/////////// Clear current denoised window
		clear(rwindow,0.0,winsize);
		clear(gwindow,0.0,winsize);
		clear(bwindow,0.0,winsize);*/

		for(int j=jmin; j <= jmax; j++)
			for(int i=imin ; i <= imax; i++) {
		
				float dif = 0.0f;


				if (kernel && nwin0==nwin){

					dif = weighted_l2_distance(ired,ired,x,y,i,j,width,kernel,nwin);
					dif += weighted_l2_distance(igreen,igreen,x,y,i,j,width,kernel,nwin);
					dif += weighted_l2_distance(iblue,iblue,x,y,i,j,width,kernel,nwin);
		

				} else
				{

					dif=l2_distance(igreen,igreen,x,y,i,j,nwin0,width);
					dif+=l2_distance(ired,ired,x,y,i,j,nwin0,width);
					dif+=l2_distance(iblue,iblue,x,y,i,j,nwin0,width);
					
				}

				dif = dif / 3.0;

				/////// If difference less than threshold average window (i,j)
				if (dif < t2)	
				{					

					float weight;

					///// weight => 1 		if dist < t1
					////		(t2-x)/(t2-t1)	if t1<dist<t2
					////		0		if t>t2
				
					//// weight function is continuous but not derivable at t1 and t2
					 
					if (dif < t1) weight = 1.0f;
					else weight = (t2 - dif) / t2mt1;	// t2mt1 = t2 - t1 (declared outside fors)


					///// We fill window of size 2*nwin0+1 x 2*nwin0+1
					for(int r=-nwin0; r<=nwin0; r++)
						for(int s=-nwin0; s <= nwin0; s++)
					{ 

						int index = (y+s)*width + x+r;	//// However the complete rwindow/gwindow/bwindow have size depending on nwin
						int l=(j+s)*width+i+r;

						ored[index] +=   weight * ired[l];
						ogreen[index] +=   weight * igreen[l];
						oblue[index] +=   weight * iblue[l];	

						counter_mask[index] += weight;
					}		

	
				}


			}
      
		}

			
		for(int i=0; i < wxh;i++)
			if (counter_mask[i]>0.0) 
			{
					ored[i] /= counter_mask[i];
					ogreen[i] /= counter_mask[i];
					oblue[i] /= counter_mask[i];
	
			} else  //// This should never happen, since at least for each pixel itself has distance 0
			{
				ored[i] = ired[i];
				ogreen[i] = igreen[i];
				oblue[i] = iblue[i];
			}

}
