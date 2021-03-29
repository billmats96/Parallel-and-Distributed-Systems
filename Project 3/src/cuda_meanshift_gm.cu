/*
 * cuda_meanshift_gm.c
 *
 *  Created on: Jan 23, 2018
 *      Author: Matsoukas Vasileios
 *
 * Implementation using global memory instead of shared.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>

__device__ double kernel_fun(double *x, double *y,int i, int j, double h,int numOfCoordinates);
__device__ double norm(double *m,int totalPoints,int numOfCoordinates);
void test(double *y,int length);
__global__ void perform_meanshift(double *x, double *y,double *y_new, double*m,double *gdata_num,double*gdata_den,double h,int totalPoints,int numOfCoordinates,int MAXROWS,int threads,int *k);
void setMAXROWS_THREADS(); //this function balances the distribution of shared memory according to the dimension of input set
void setMAXROWS_THREADS2();//redistribution if set has less than 1024 points

struct timeval startwtime, endwtime;
double seq_time;
double *matlabY;
int MAXROWS,THREADS,numOfCoordinates,totalPoints;

int main(int argc, char **argv)
{
	FILE *fp;
	int sz,k=0,length,threadsPerBlock,blocksPerGrid;  //host variables
	double *x, *y, *y_new,*m,*gdata_num,*gdata_den;  //host variables

	double *d_x,*d_y,*d_y_new,*d_m; //device copies
	int *d_k; //device copies
	int testMode=0; //Set testMode 0 to use your own input. Set testMode to 1 if you want to run tests using random inputs with matlab script.

	if (argc!=6)
	{
			printf("wrong number of inputs..");
			exit(0);
	}

	//Define h, number of iterations, dimensions of set
	double h=atof(argv[2]);
	int iterations=atoi(argv[3]);
	numOfCoordinates=atoi(argv[5]);

    setMAXROWS_THREADS();

	//Open data file

	fp=fopen(argv[1],"rb");

	if (!fp)
	{
		printf("Unable to open file!");
		return 1;
	}

	//Extract file total size and divide by sizeof(double) to find total data length
	//Then divide with number of coordinates of each point, to find out file's total points.
	fseek(fp, 0, SEEK_END);
	sz = ftell(fp);
	printf("Size of data.bin in bytes=%d\n",sz);
	length= sz/sizeof(double); // take size in groups of double
	printf("Length=%d\n",length);
	x=(double *)malloc(length*sizeof(double)); //Dynamically allocate space in memory
	printf("Size of double = %d bytes\n",sizeof(double));

	if(x==NULL)
	{
		printf("Cannot allocate memory..\n");
		exit(1);
	}

	rewind(fp);
	fread(x,sizeof(double),length,fp);
	fclose(fp);

	//Print your input data
//	printf("\nThat's the data inserted:");
//	for (int i=0; i<length ; i++)  //change to -->length
//	{
//		if(i%numOfCoordinates==0) printf("\n");
//		printf("%f ",x[i]);
//	}
//	printf("\n\n");

	totalPoints=length/numOfCoordinates;

	y=(double *)malloc(length*sizeof(double));
	for(int i=0 ; i<length; i++)
	{
		y[i]=x[i];
	}
//	printf("Initial y\n");
//	for (int i=0; i<length; i++)  //change to -->length
//	{
//		if(i%numOfCoordinates==0) printf("\n");
//		printf("%f ",y[i]);
//	}
	y_new=(double *)malloc(length*sizeof(double));
	m=(double *)malloc(length*sizeof(double));
	printf("\n\n");

	//Intialize y_new to zero and m to infinity.
	for (int i=0; i<length ; i++)
	{
	    y_new[i]=0;
	    m[i]=INFINITY;
	}

	if (totalPoints<1024)
	{
		setMAXROWS_THREADS2();
	}
	//Allocate space for device copies
	cudaMalloc(&d_x, length*sizeof(double));
	cudaMalloc(&d_y, length*sizeof(double));
	cudaMalloc(&d_y_new, length*sizeof(double));
	cudaMalloc(&d_m, length*sizeof(double));
	cudaMalloc(&d_k,sizeof(int));
	//Allocate same portion of global memory as total shared used in the other implementation
	cudaMalloc(&gdata_num,MAXROWS*totalPoints*THREADS*sizeof(double));
	cudaMalloc(&gdata_den,totalPoints*THREADS*sizeof(double));

	//Copy inputs to device
	cudaMemcpy(d_x, x, length*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_y, y, length*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_y_new, y_new, length*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_m, m, length*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_k, &k, sizeof(int), cudaMemcpyHostToDevice);

	//Define Blocks and Threads/Block to use
	threadsPerBlock = THREADS;
	blocksPerGrid   = totalPoints;

	//Begin chronometer
	gettimeofday (&startwtime, NULL);

	//let user define constant iterations. Perform meanshift.
	for(int t=0;t<iterations;t++)
	{
		perform_meanshift<<<blocksPerGrid,threadsPerBlock>>>(d_x,d_y,d_y_new,d_m,gdata_num,gdata_den,h,totalPoints,numOfCoordinates,MAXROWS,THREADS,d_k);
	}

	//Stop chronometer
	cudaDeviceSynchronize();
	gettimeofday (&endwtime, NULL);

	//Copy results back to host
	cudaMemcpy(&k, d_k, sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(y, d_y, length*sizeof(double), cudaMemcpyDeviceToHost);

	//Print output data
	//printf("\nThe y_final i calculated:\n");
	//for(int i=(length); i<length; i++)
	//{
	//  if(i%numOfCoordinates==0) printf("\nPoint %d: ",i/numOfCoordinates);
	//   printf("%f ",y[i]);
	//}
	printf("\n");

	seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6+ endwtime.tv_sec - startwtime.tv_sec);
	printf("Total time needed to execute mean shift using global memory: %f secs\n", seq_time);
	printf("Needed %d iterations\n",k);

	//Write results to a .bin file.
	fp=fopen(argv[4],"wb");
	fwrite(y,length,sizeof(double),fp);
	fclose(fp);

	//Test results if test mode is true. Use it only with matlab script
	if(testMode)
	{
		matlabY=(double *)(malloc(length*sizeof(double)));

		fp=fopen("outputmatlab.bin","rb");
		fread(matlabY,sizeof(double),length,fp);
		fclose(fp);
		test(y,length);
	}


	//Free allocated memory
	free(x);
	free(y);
	free(y_new);
	free(m);
	free(matlabY);

	cudaFree(d_x);
	cudaFree(d_y);
	cudaFree(d_y_new);
	cudaFree(d_m);
	cudaFree(d_k);
	cudaFree(gdata_num);
	cudaFree(gdata_den);

	return 0;
}


//Gaussian Kernel Function.Take note that we calculate distance between points without
//taking the square root of it.
__device__ double kernel_fun(double *x, double *y,int i, int j, double h,int n)
{
	double distance=0;
	for(int l=0; l<n;l++)
	{
		distance+=(y[n*i+l]-x[n*j+l])*(y[n*i+l]-x[n*j+l]);
	}
	if (distance<=h*h)
	{
		return exp(-distance/(2*h*h));
	}
	else
	{
		return 0;
	}
}

//Frobenius norm calculation of mean_shift  "m"  matrix
__device__ double norm(double *m,int totalPoints,int n)
{
	double norm_value=0;
	for(int i=0; i<totalPoints; i++)
	{
		for(int l=0; l<n;l++)
			{
				norm_value+=m[n*i+l]*m[n*i+l];
			}
	}
	return sqrt(norm_value);
}

void test(double *y,int length)
{
	int flag=0;
	int counter=0;
	for(int i=0; i<length; i++)
	{
		if(fabs(y[i]-matlabY[i])>=0.2)
		{
			counter++;
		}
	}
	if((double)counter/(double)length >0.02)
	{
		printf("%d coordinates that differ more than 0.2\n",counter);
		printf("Considerable error\n");
		flag=1;
	}
	if(flag)
	{
		printf("Initial Test Failed!\n");
		printf("Check norm value of comparison matrix C\n");
	}
	else
	{
		printf("Test Passed\n");
	}
}

__global__ void perform_meanshift(double *x, double *y,double *y_new, double*m,double *gdata_num,double*gdata_den,double h,int totalPoints,int numOfCoordinates,int MAXROWS,int threads,int *k)
{
	double den=0; //holds the value of denominator of each yi
	double temp_dist=0; //store the result of kernel_fun() to avoid multiple function calls
	int n=numOfCoordinates; //for convenience
	int i=blockIdx.x;
	int j =threadIdx.x;
	int z=blockDim.x*blockIdx.x + threadIdx.x;
	int temp=0;

	  if(z==0)*k=*k+1; //count iterations in thread thread 0, block 0
		if(i<totalPoints)
		{
			if(j<totalPoints)
			{
				//l is the current dimension. Calculate MAXROWS dimensions in each loop
				for(int l=0 ; l<n ; l+=MAXROWS)
				{
					//temp is used to define the right number of dimensions to process.
					if((n-l)<MAXROWS)
					{
						temp=n-l;
					}
					else
					{
						temp=MAXROWS;
					}
					//Parallelize the sum of numerator and denominator.each thread writes in its own position in global memory
					//Set initial content to zero
					//#pragma unroll
					for(int c=0; c<temp; c++)
					{
						gdata_num[j+threads*(c+MAXROWS*i)]=0;
					}
					if (l==0) gdata_den[j+threads*i]=0;

					//Each thread takes one point. In case points>1024 (or max num of threads) then each thread takes up some more points
					for(int q=j ;q<totalPoints;q+=threads)
					{
						temp_dist=kernel_fun(x,y,i,q,h,n);
					//#pragma unroll
						for(int c=0; c<temp; c++)
						{
							gdata_num[j+threads*(c+MAXROWS*i)]+=temp_dist*x[n*q+l+c]; //all threads do in parallel the multiplication. for every coordinate (l+c)
						}
						if(l==0) gdata_den[j+threads*i]+=temp_dist; //denonimator is once calculated for each yi
					}
					__syncthreads();
					//Now do the reduction to sum threads' values and have numerator's value
					//Threads have to be power of 2 for that reason.The program takes care to assign write number of threads.
					for (unsigned int s=blockDim.x/2; s>0; s>>=1)
					{
					    if (j<s)
						{
					    	for(int c=0; c<temp; c++)
					    	{
					    		gdata_num[j+threads*(c+MAXROWS*i)]+=gdata_num[j+threads*(c+MAXROWS*i)+s];
					    	}
					    	if(l==0) gdata_den[j+threads*i]+=gdata_den[j+threads*i+s];
						}
					    __syncthreads();
					}
					if(l==0) den=gdata_den[0+threads*i];
					if(j==0)
					{
						for(int c=0; c<temp; c++)
						{
							y_new[n*i+l+c]=gdata_num[0+threads*(c+MAXROWS*i)];
						}

					}
					__syncthreads();
				}
			}
			__syncthreads();
			if (j == 0 && den!=0)
			{
				for(int l=0 ; l<n ; l++)
				{
					y_new[n*i+l]=y_new[n*i+l]/den;
					//m[n*i+l]=y_new[n*i+l]-y[n*i+l];
				}
				for(int l=0 ; l<n; l++)
			    {
					y[n*i+l]=y_new[n*i+l];
					y_new[n*i+l]=0;
				}
			}
			__syncthreads();
			den=0; //reset den after each point iteration
		}

}
void setMAXROWS_THREADS()
{
	if (numOfCoordinates<=5)
	{
		MAXROWS=5;
		THREADS=1024;
	}
	else if(numOfCoordinates>5 && numOfCoordinates<=11)
	{
		MAXROWS=11;
		THREADS=512;
	}
	else if(numOfCoordinates>11 && numOfCoordinates<=23)
	{
		MAXROWS=23;
		THREADS=256;
	}
	else if(numOfCoordinates>23 && numOfCoordinates<=47)
	{
		MAXROWS=47;
		THREADS=128;
	}
	else
	{
		MAXROWS=95;
		THREADS=64;
	}
}

void setMAXROWS_THREADS2()
{
	int n=(int)ceil(THREADS/totalPoints);
	int i=1;
	while (n>pow(2,i))
	{
		i++;
	}
	n=pow(2,i);
	THREADS=THREADS/n;
	MAXROWS=n*MAXROWS+n-1;
}


