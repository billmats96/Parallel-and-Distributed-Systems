/*
 * Non Local Means - Parallel Implementation
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
#include <string.h>

#define MAX 512

struct timeval startwtime, endwtime;
float seq_time=0;
float* image_pad(float *im,int m, int n,int patchsize_x,int patchsize_y);
float* gauss_patch(float* gaussian,int offsetx, int offsety, float patchSigma);


//------------Kernel functions------------//


// finding each pixel's neighbors

__global__ void neighbs(float* im, float* neighbors, float* gauss, int offsety, int offsetx, int patch_size, int newpatch_size, int s_x, int s_y, int patchsize_y){
	int id = threadIdx.x+blockIdx.x*blockDim.x;
	int m = s_y - 2*offsety;
	int n = s_x - 2*offsetx;
	if (id<m*n){
		int kx = blockIdx.y-offsetx;
		int ky = blockIdx.z-offsety;
		int neighbor_id = blockIdx.y*patchsize_y+blockIdx.z;
		neighbors[id*newpatch_size+neighbor_id] = im[(offsety+ky+id%m)*s_x+kx+offsetx+(id/m)]*gauss[neighbor_id];

	}
}




__global__ void affinity(float *image, float *neighbors, int pixels, int patchsize, float filtSigma, int offsetx,int offsety,int patchsize_y, float *row){
	__shared__ float Ys[16][16];
	__shared__ float Xs[16][16];

	int bx = blockIdx.x, by = blockIdx.y;
	int tx = threadIdx.x, ty = threadIdx.y;

	int yBegin = by * 16 * patchsize;
	int xBegin = bx * 16 * patchsize;
	int yEnd = yBegin + patchsize - 1, y, x, k;
	float tmp, c = 0 , s = 0;

	int pixel_x = bx*16 + tx;
	int pixel_y = by*16 + ty;

	for(y=yBegin,x=xBegin; y<=yEnd; y+=16,x+=16){

		if (pixel_y < pixels) Ys[ty][tx] = neighbors[y + ty*patchsize + tx];
		else Ys[ty][tx] = 0;
		if (bx*16+ty < pixels) Xs[tx][ty] = neighbors[x + ty*patchsize + tx];
		else Xs[tx][ty] = 0;

		__syncthreads();

		for(k=0;k<16;k++){
			tmp = Ys[ty][k] - Xs[k][tx];
			s += tmp*tmp;
		}
		__syncthreads();
	}

	if (pixel_y < pixels && pixel_x < pixels){
		Xs[ty][tx] = exp(-s/(filtSigma*filtSigma));
		Ys[ty][tx] = Xs[ty][tx]*neighbors[pixel_x*patchsize+patchsize_y*offsetx+offsety];
	}
	else {
		Ys[ty][tx] = 0;
		Xs[ty][tx] = 0;
	}
	__syncthreads();

	if (pixel_y< pixels && tx==0){
		s=0;
		for (k=0; k<16; k++){
			c+=Ys[ty][k];
			s+=Xs[ty][k];
		}
		atomicAdd(&image[pixel_y], c); atomicAdd(&row[pixel_y],s);
	}


}


__global__ void newimage(float* im, float* row, int k, int n){
	int id = blockIdx.x*n + blockIdx.y;
	if (id < k) im[id] = im[id]/row[id];
}



int main(int argc, char **argv){

	if (argc!=5){
		printf("Wrong number of inputs.\n");
		printf("Usage: %s binary_input binary_output image_rows image_columns \n",argv[0]);
		exit(0);
	}


	FILE *fp;					//file pointer to the binary image file
	int m = atoi(argv[3]); 				//image rows
	int n = atoi(argv[4]); 				//image columns


	int patchsize_x = 5, patchsize_y = 5;
	int patch_size = patchsize_x*patchsize_y;
	float patchSigma = 1.6667;
	float filtSigma = 0.1;

	int size_y = m + (patchsize_y-1);
	int offsety = (patchsize_y-1)/2;
	int size_x = n + (patchsize_x-1);
	int offsetx = (patchsize_x-1)/2;

	float *im, *new_im, *gaussian;
	float *d_im, *d_new_im, *d_gaussian, *d_neighbors, *d_row;
	int blocksx, blocksy, blocksz;

	//Allocate (1-D) memory for input image and filtered image

	im = (float *)malloc(size_y*size_x*sizeof( float ));
	new_im = (float *)malloc(m*n*sizeof(float));

	//Open data file

	fp=fopen(argv[1],"rb");

	if (!fp){
		printf("Unable to open file!");
		return 1;
	}

	//Store input image data

	for (int i=offsety; i<size_y-offsety; i++){
    		fread(im+i*size_x+offsetx,sizeof(float),n,fp);
	}

	fclose(fp);

	//Pad image's borders symmetrically

	im = image_pad(im, m, n, patchsize_x, patchsize_y);

	//Create gaussian patch

	int k = m*n;
	gaussian = (float *)malloc(patch_size*sizeof( float));
	gaussian = gauss_patch(gaussian, offsetx, offsety, patchSigma);

	// Set new patchsize

	int newpatch_size = patch_size;
	if (patch_size%16!=0) newpatch_size = (patch_size/16 + 1)*16;


	//Allocate space for device copies

	cudaMalloc((void **)&d_im, size_x*size_y*sizeof(float));
	cudaMalloc((void **)&d_new_im, m*n*sizeof(float));
	cudaMalloc((void **)&d_neighbors, (m*n)*newpatch_size*sizeof(float));
	cudaMalloc((void **)&d_gaussian, patch_size*sizeof(float));
	cudaMalloc((void **)&d_row, (m*n)*sizeof(float));


	//Copy inputs to device

	cudaMemcpy(d_im, im, size_x*size_y*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_gaussian, gaussian, patch_size*sizeof(float), cudaMemcpyHostToDevice);

	free(gaussian); free(im);

	// Cuda Memset

	cudaMemset(d_neighbors, 0, (m*n)*newpatch_size*sizeof(float));
	cudaMemset(d_row, 0, (m*n)*sizeof(float));
	cudaMemset(d_new_im, 0, (m*n)*sizeof(float));

	//Finding neighbors, number_of_neighbors % 16 = 0

	blocksx = (m*n)/MAX;
	if ((m*n)%MAX!=0)
		blocksx++;

	blocksy = patchsize_x;
	blocksz = patchsize_y;

	//Begin Chronometer
	gettimeofday (&startwtime, NULL);

	neighbs<<<dim3(blocksx,blocksy,blocksz),MAX>>>(d_im,d_neighbors, d_gaussian, offsety, offsetx, patch_size, newpatch_size, size_x, size_y, patchsize_y);

	//Stop chronometer
	cudaDeviceSynchronize();
	gettimeofday (&endwtime, NULL);

	printf("\n");
	seq_time += (float)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6+ endwtime.tv_sec - startwtime.tv_sec);
	//printf("Total time needed to find neighbors using gpu : %f secs\n", seq_time);



	//Finding the affinity matrix, k%16 = 0

	blocksx = k/16;

	if (k%16!=0) blocksx = (blocksx/16+1)*16;
	blocksy = blocksx;

	//Begin chronometer
 	gettimeofday (&startwtime, NULL);

	affinity<<<dim3(blocksx, blocksy, 1),dim3(16, 16, 1)>>>(d_new_im, d_neighbors, k, newpatch_size, filtSigma, offsetx, offsety, patchsize_y, d_row);

	newimage<<<dim3(m, n), 1>>>(d_new_im, d_row, k, n);

    	//Stop chronometer
	cudaDeviceSynchronize();
	gettimeofday (&endwtime, NULL);


	printf("\n");
	seq_time += (float)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6+ endwtime.tv_sec - startwtime.tv_sec);
	//printf("Total time needed to create affinity matrix D using gpu : %f secs\n", seq_time);
	printf("Total time needed using nlmeans: %f secs\n", seq_time);

    	cudaMemcpy(new_im,d_new_im, m*n*sizeof(float), cudaMemcpyDeviceToHost );

//    	printf("Thats the filtered image:\n");
//
//
//		for (int j=0; j<m; j++){
//			for (int i=0; i<n ; i++)
//	       			printf( "%f ",new_im[i*m+j]);
//
//	  	        printf("\n\n");
//	    }

		//Save output to a binary file
		fp=fopen(argv[2],"w");
		for (int j=0; j<m; j++){
			for (int i=0; i<n ; i++)
				fwrite(&new_im[i*m+j],sizeof(float),1,fp);
	    }
		fclose(fp);

    	//Free allocated memory

	cudaFree(d_im); cudaFree(d_gaussian); cudaFree(d_neighbors); cudaFree(d_row); cudaFree(d_new_im);

    	free(new_im);

    	return 0;

}






//------------Host functions------------//

float* image_pad(float *im,int m, int n,int patchsize_x,int patchsize_y){
	int size_y = m + (patchsize_y-1);
	int offsety = (patchsize_y-1)/2;
	int size_x = n + (patchsize_x-1);
	int offsetx = (patchsize_x-1)/2;

	//Padding the data for image's borders

	for (int i=0; i<offsety; i++){
		for (int j=offsetx ; j<size_x - offsetx; j++)
			im[i*size_x+j] = im[(2*offsety-1-i)*size_x+j];
	}

	for (int i=0; i<offsety; i++){
		for (int j=offsetx ; j<size_x - offsetx; j++)
			im[(size_y-offsety+i)*size_x+j] = im[(size_y-offsety-i-1)*size_x+j];
	}

	for (int i=0; i<size_y; i++){
		for (int j=0 ; j<offsetx; j++)
			im[i*size_x+j] = im[i*size_x+2*offsetx-1-j];
	}

	for (int i=0; i<size_y; i++){
		for (int j=0 ; j<offsetx; j++)
			im[i*size_x+size_x-offsetx+j] = im[i*size_x+size_x-offsetx-1-j];
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

