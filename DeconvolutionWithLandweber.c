#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "Function.h"


#define NAME_IMG_IN  "photograph"

#define NAME_IMG_OUT1 "photograph_original"
#define NAME_IMG_OUT2 "photograph_degraded_withoutNoise" 
#define NAME_IMG_OUT3 "photograph_restored_withoutNoise"  
#define NAME_IMG_OUT4 "photograph_degraded_withNoise" 
#define NAME_IMG_OUT5 "photograph_restored_withNoise"  
#define alpha 1

int main(int argc,char** argv)
{
	int i,j,k;
  	int nb_iterations;
  	int length,width;
  	float var, isnr;
	float sum1,sum2;
  	int size_filter; 

  	float** image; /* original image */
  	float** g;  /* degraded image */
  	float** f;  /* restored image */
        float** h;  /* blur function*/
   	float** temp;

  	printf("Input size of low-pass filter : ");
  	scanf("%d",&size_filter);

  	printf("\nInput the number of iteration for LANDWEBER: ");
  	scanf("%d",&nb_iterations);
		
	image=LoadImagePgm(NAME_IMG_IN, &length, &width);
	g=fmatrix_allocate_2d(length,width);
	f=fmatrix_allocate_2d(length,width);
	h=fmatrix_allocate_2d(length,width);
	temp=fmatrix_allocate_2d(length,width);

	for(i=0;i<length;i++){
		for(j=0;j<width;j++){
			g[i][j]=0.0;
			f[i][j]=0.0;
			h[i][j]=0.0;	
			temp[i][j]=0.0;
		}
	}

  	/* add blur to original image : g = image + blur */
	flou(h,size_filter,length,width);
	G(g,image,h,length,width);

 	/*****************************************/
	/* restore image g with LANDWEBER.       */
	/* and display ISNR for every iteration  */
	/*****************************************/
	copy(f,g,length,width);
	for(k=1;k<=nb_iterations;k++){
		sum1=0.0;
		sum2=0.0;
		for(i=0;i<length;i++){
			for(j=0;j<width;j++){
				sum1=sum1+SQUARE(image[i][j]-g[i][j]);
				sum2=sum2+SQUARE(image[i][j]-f[i][j]);
			}
		}

		printf("\n %d\t %f\t",k,10*log10(sum1/sum2));
		update(f,temp,h,g,alpha,size_filter,length,width);
		
	}
	
	Recal(f,length,width);

  	/*Save images */
	SaveImagePgm(NAME_IMG_OUT1,image,length,width);
	system("display photograph_original.pgm&");
        SaveImagePgm(NAME_IMG_OUT2,g,length,width);
	system("display photograph_degraded_withoutNoise.pgm&");
	SaveImagePgm(NAME_IMG_OUT3,f,length,width);
	system("display photograph_restored_withoutNoise.pgm&");


  	printf("Input variance of gaussian noise: ");
  	scanf("%f",&var);
	
  	/* Add noise and blur to image:  image = image + blur + noise (add_gaussian_noise) */
	flou(h,size_filter,length,width);
	G(g,image,h,length,width);
	add_gaussian_noise(g,length,width,var);
	
  	/************************************************/
	/* restore image with LANDWEBER.                */
	/* and display ISNR for every iteration         */
	/************************************************/
	copy(f,g,length,width);
	for(k=1;k<=nb_iterations;k++){
		sum1=0.0;
		sum2=0.0;
		for(i=0;i<length;i++){
			for(j=0;j<width;j++){
				sum1=sum1+SQUARE(image[i][j]-g[i][j]);
				sum2=sum2+SQUARE(image[i][j]-f[i][j]);
			}
		}
		printf("\n %d\t %f\t",k,10*log10(sum1/sum2));
		update(f,temp,h,g,alpha,size_filter,length,width);
		
	}
	
	Recal(f,length,width);	

  	/*Save images */
	SaveImagePgm(NAME_IMG_OUT4,g,length,width);
	system("display photograph_degraded_withNoise.pgm&");
	SaveImagePgm(NAME_IMG_OUT5,f,length,width);
	system("display photograph_restored_withNoise.pgm&");

  	/*free memory*/
	free_fmatrix_2d(image);  
	free_fmatrix_2d(f);
	free_fmatrix_2d(g);
	free_fmatrix_2d(h);
	free_fmatrix_2d(temp);
  	
  	printf("\n Ending \n\n\n");
  	return 0; 	 
}
