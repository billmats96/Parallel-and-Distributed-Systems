/*
 * Adaptive Non Local Means - Serial Implementation
 *
 * Authors: Athanasiadis Christos athanasc@ece.auth.gr AEM 8416
 *     	    Matsoukas Vasileios   vmatsouk@ece.auth.gr AEM 8743
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/time.h>
#include <math.h>


float** image_pad(float** im,int m, int n,int patchsize_x, int patchsize_y);
float* gauss_patch(float *gaussian, int offsetx,int offsety,float patchSigma);
float** neighbs(float** neighbors, float** im,float* gaussian,int m, int n, int offsetx, int offsety);
float* filter_image(float* newim, float** neighbors,int *regions, int m, int n,int patch_size,int patchsize_x,int offsetx, int offsety,int regionscount,float std);

struct timeval startwtime, endwtime;
double seq_time=0;


int main(int argc, char **argv){

	if (argc!=5){
		printf("Wrong number of inputs.\n");
		printf("Usage: %s binary_input binary_output image_rows image_columns \n",argv[0]);
		exit(0);
	}

	FILE *fp;
	int m = atoi(argv[3]);
	int n = atoi(argv[4]);


	int patchsize_x = 5, patchsize_y = 5;
	int patch_size = patchsize_x*patchsize_y;
	float patchSigma = 1.6667; float filtSigma=0.1;
	int nLevel = 6;

	int size_y = m + (patchsize_y-1);
	int offsety = (patchsize_y-1)/2;
	int size_x = n + (patchsize_x-1);
	int offsetx = (patchsize_x-1)/2;

	float *gaussian,**neighbors,**im,*newim,*mean,*std;
	int **L, *count,**regions,*regionscount,reg,counter;

	im = ( float **)malloc(size_y*sizeof( float *));
	for (int i=0; i<size_y; i++)
		im[i] = ( float *)malloc(size_x*sizeof( float));


	//Open data file
	fp=fopen(argv[1],"rb");

	if (!fp){
		printf("Unable to open file!");
		return 1;
	}


	//Store input data
    for (int i=offsety; i<size_y-offsety; i++){
    	fread(*(im+i)+offsetx,sizeof(float),n,fp);
    }

	fclose(fp);



	//Distinct Levels
	L = ( int **)malloc(m*sizeof( int *));
	for (int i=0; i<m; i++)
		L[i] = ( int *)malloc(n*sizeof( int ));

	mean = (float *)calloc(nLevel,sizeof( float ));
	std = (float *)calloc(nLevel,sizeof( float ));
	count = (int *)calloc(nLevel,sizeof(int));

	for (int j=0; j<n; j++){
		for (int i=0; i<m; i++){
			L[i][j] = round((nLevel-1)*im[i+offsety][j+offsetx]);
			mean[L[i][j]]+=im[i+offsety][j+offsetx];
			count[L[i][j]]++;
		}
	}

//	//Print Level values
//	printf("The level of each pixel: \n");
//	for (int i=0; i<m; i++){
//		for (int j=0; j<n; j++){
//			printf("%d ",L[i][j]);
//		}
//		printf("\n");
//	}
//	printf("\n\n");


	for (int i=0; i<nLevel; i++) mean[i]/=count[i];


	//Create regions with pixels of same level
	regions = ( int **)malloc(nLevel*sizeof( int *));
	for (int i=0; i<nLevel; i++)
		regions[i] = ( int *)malloc(count[i]*sizeof(int));

	regionscount = (int *)calloc(nLevel,sizeof(int));


    for (int j=0; j<n; j++){
		for (int i=0; i<m; i++){
			reg=L[i][j];
			regions[reg][regionscount[reg]++] = j*m+i;
			std[L[i][j]]+=(im[i+offsety][j+offsetx]-mean[L[i][j]])*(im[i+offsety][j+offsetx]-mean[L[i][j]]);
		}
    }

    for (int i=0; i<nLevel; i++) std[i]=sqrt(std[i]/(count[i]-1));

//	//Print std value for each region
//  printf("The standard deviation value of each region: \n");
//	for (int i=0; i<nLevel; i++){
//		printf("%f  ", std[i]);
//	}
//	printf("\n\n");
//
//	//Print the pixel IDs of each region
//	printf("The pixels of each region: \n");
//	for (int i=0; i<nLevel; i++){
//		for(int j=0; j<count[i]; j++){
//			printf("%d ",regions[i][j]);
//		}
//		printf("\n\n");
//	}
//	printf("\n\n");


	//Pad image's borders symmetrically
	im=image_pad(im,m,n,patchsize_x,patchsize_y);


	//Create gaussian patch
	gaussian = ( float *)malloc(patch_size*sizeof( float));
	gaussian = gauss_patch(gaussian,offsetx,offsety,patchSigma);


	//Finding neighbors (im2col matlab)
	neighbors = ( float **)malloc(m*n*sizeof( float *));
	for (int i=0; i<m*n; i++)
		neighbors[i] = ( float *)malloc(patch_size*sizeof( float));

	gettimeofday (&startwtime, NULL);
	neighbors = neighbs(neighbors,im,gaussian,m,n,offsetx,offsety);
    gettimeofday (&endwtime, NULL);

	printf("\n");
	seq_time += (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6+ endwtime.tv_sec - startwtime.tv_sec);
	//printf("Total time needed to find neighbors using cpu : %f secs\n", seq_time);


	//Image Filtering
	newim = ( float *)calloc(m*n,sizeof( float *));

	gettimeofday (&startwtime, NULL);
	for(int i=0; i<nLevel; i++){

		newim=filter_image(newim,neighbors,regions[i],m,n,patch_size,patchsize_x,offsetx,offsety,count[i],std[i]);

	}
	gettimeofday (&endwtime, NULL);

	printf("\n");
	seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6+ endwtime.tv_sec - startwtime.tv_sec);
	//printf("Total time needed to filter image using cpu : %f secs\n\n", seq_time);
	printf("Total time needed using serial adaptive nlmeans: %f secs\n", seq_time);

//	//Print filtered image
//	printf("That's the filtered image :\n");
//	counter=0;
//	for(int i=0; i<m*n; i++){
//
//		 printf( "%f ", newim[i/n + (i%n)*m]);
//		 counter++;
//		 if (counter==n){
//			 printf("\n");
//			 counter=0;
//		 }
//
//	}

	//Save output to a binary file
	fp=fopen(argv[2],"w");
	for(int i=0; i<m*n; i++){
		 fwrite(&newim[i/n + (i%n)*m],sizeof(float),1,fp);
	}
	fclose(fp);

	//Free allocated memory
	for (int i=0; i<size_y; i++){
		free(im[i]);
	}
	for (int i=0; i<m*n; i++){
		free(neighbors[i]);
	}
	for (int i=0; i<m; i++){
		free(L[i]);
	}
	for (int i=0; i<nLevel; i++){
		free(regions[i]);
	}

	free(im);
	free(gaussian);
	free(neighbors);
	free(newim);
	free(mean);
	free(std);
	free(count);
	free(L);
	free(regionscount);
	free(regions);

	return 0;


}

float** image_pad(float** im,int m, int n,int patchsize_x, int patchsize_y)
{
	int size_y = m + (patchsize_y-1);
	int offsety = (patchsize_y-1)/2;
	int size_x = n + (patchsize_x-1);
	int offsetx = (patchsize_x-1)/2;


	for (int i=0; i<offsety; i++){
		for (int j=offsetx ; j<size_x - offsetx; j++)
			im[i][j] = im[2*offsety-1-i][j];
	}
	for (int i=0; i<offsety; i++){
		for (int j=offsetx ; j<size_x - offsetx; j++)
			im[size_y-offsety+i][j] = im[size_y-offsety-i-1][j];
	}

	for (int i=0; i<size_y; i++){
		for (int j=0 ; j<offsetx; j++)
			im[i][j] = im[i][2*offsetx-1-j];
	}
	for (int i=0; i<size_y; i++){
		for (int j=0 ; j<offsetx; j++)
			im[i][size_x-offsetx+j] = im[i][size_x-offsetx-1-j];
	}

	return im;
}

float* gauss_patch(float *gaussian, int offsetx,int offsety,float patchSigma){
	int u=0;

	for (int kx=-offsetx; kx<offsetx+1; kx++){
		for (int ky = -offsety; ky<offsety+1; ky++){
			gaussian[u] = exp(-(ky*ky+kx*kx)/(2*patchSigma*patchSigma));
			u++;
		}
	}

	return gaussian;
}

float** neighbs(float** neighbors, float** im,float* gaussian,int m, int n, int offsetx, int offsety)
{
	int id=0;
	int neighbor_id=0;

	for (int j=0; j<n; j++){
		for (int i=0; i<m; i++){
			for (int kx=-offsetx; kx<offsetx+1; kx++){
				for (int ky = -offsety; ky<offsety+1; ky++){
					neighbors[id][neighbor_id] = im[i+offsety+ky][j+offsetx+kx]*gaussian[neighbor_id];
					neighbor_id++;
				}
			}
			neighbor_id=0;
			id++;
		}
	}


	return neighbors;
}

float* filter_image(float* newim, float** neighbors,int *regions, int m, int n,int patch_size,int patchsize_x,int offsetx, int offsety,int regionscount,float std)
{
	float filtVar = std*std;
	float dist;
	int i,j;

	float *sumofrow=(float *)calloc(m*n,sizeof(float));

	for (int q=0; q<regionscount; q++){
		i=regions[q];
		for (int v=q; v<regionscount; v++){
			j=regions[v];
			dist=0.0;

			for (int d=0; d<patch_size; d++)
				dist+=(neighbors[i][d]-neighbors[j][d])*(neighbors[i][d]-neighbors[j][d]);

			dist = exp(-dist/filtVar);
			sumofrow[i]+=dist;
			newim[i]+=dist*neighbors[j][patchsize_x*offsety+offsetx];
			if(i!=j){
				newim[j]+=dist*neighbors[i][patchsize_x*offsety+offsetx];
				sumofrow[j]+=dist;
			}

		}
	}

	for (int q=0; q<regionscount; q++){
		i=regions[q];
		newim[i]=newim[i]/sumofrow[i];
	}

	return newim;
}




