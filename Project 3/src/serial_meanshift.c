/*
 * serial_meanshift.c
 *
 *  Created on: Jan 18, 2018
 *      Author: linuxbill
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>

double kernel_fun(double *x, double *y,int i, int j, double h,int numOfCoordinates);
double norm(double *m,int totalPoints,int numOfCoordinates);
void test(double *y,int length);
struct timeval startwtime, endwtime;
double seq_time;
double *matlabY;

int main(int argc, char **argv)
{
	FILE *fp;
	int sz,k=0,totalPoints,length; //value of numOfCoordinates according to the coordinates of dataset to be used.
	double *x, *y, *y_new,*m,epsilon=0.0001,h;
	double den;
	int numOfCoordinates;
	int testMode=0; //Set testMode 0 to use your own input. Set testMode to 1 if you want to run tests using random inputs with matlab script.


	//Open data file
	fp=fopen(argv[1],"rb");
	if (!fp)
	{
		printf("Unable to open file!");
		return 1;
	}

	h=atof(argv[2]);
	numOfCoordinates=atoi(argv[4]);

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

	//Begin chronometer
	gettimeofday (&startwtime, NULL);

	//Perform mean shift till norm(m) converges
	while(norm(m,totalPoints,numOfCoordinates)>epsilon)
	{
		int n=numOfCoordinates; //just for convenience
		k=k+1;
		for(int i=0; i<totalPoints; i++)
		{
			for(int j=0; j<totalPoints; j++)
			{
				for(int l=0 ; l<n ; l++)
				{
					y_new[n*i+l]+=kernel_fun(x,y,i,j,h,n)*x[n*j+l];
				}
			}
			for(int j=0; j<totalPoints; j++)
			{
				den+=kernel_fun(x,y,i,j,h,n);

			}
			for(int l=0 ; l<n ; l++)
			{
				y_new[n*i+l]=y_new[n*i+l]/den;
				m[n*i+l]=y_new[n*i+l]-y[n*i+l];
			}
			den=0; //reset den after each point iteration
		}
//		printf("Initial y_new\n");
//		for (int i=0; i<(length/50) ; i++)  //change to -->length
//		{
//			if(i%numOfCoordinates==0) printf("\n");
//			printf("%f ",y_new[i]);
//		}
		printf("\n");
		for(int p=0 ; p<length; p++)
	    {
			y[p]=y_new[p];
			y_new[p]=0;
		}
//		printf("Mean shift: \n");
//		for (int i=0; i<(length/50) ; i++)  //change to -->length
//		{
//			if(i%numOfCoordinates==0) printf("\n");
//			printf("%f ",m[i]);
//		}
//		sleep(1);
		printf("Iteration %d error %f",k,norm(m,totalPoints,numOfCoordinates));

	}
	//Stop chronometer
	gettimeofday (&endwtime, NULL);

	//Print output data
//	printf("\nThe y_final i calculated:\n");
//	for(int i=0; i<length; i++)
//	{
//	  if(i%numOfCoordinates==0) printf("\nPoint %d: ",i/numOfCoordinates);
//	  printf("%f ",y[i]);
//	}

	printf("\n");
	seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6+ endwtime.tv_sec - startwtime.tv_sec);
	printf("Total time needed to execute serial mean shift: %f secs\n", seq_time);
	printf("Needed %d iterations\n",k);

	//Write results to a .bin file.
	fp=fopen(argv[3],"wb");
	fwrite(y,length,sizeof(double),fp);
	fclose(fp);

	if(testMode)
	{
		matlabY=(double *)(malloc(length*sizeof(double)));

		fp=fopen("outputmatlab.bin","rb");
		fread(matlabY,sizeof(double),length,fp);
		fclose(fp);
		test(y,length);
	}
	//Test serial knn.set testMode to true to run tests with matlab script



	//Free allocated memory
	free(x);
	free(y);
	free(y_new);
	free(m);
	free(matlabY);
	return 0;
}


//Gaussian Kernel Function.
double kernel_fun(double *x, double *y,int i, int j, double h,int numOfCoordinates)
{
	double distance=0;
	int n=numOfCoordinates;
	for(int l=0; l<numOfCoordinates;l++)
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
double norm(double *m,int totalPoints,int numOfCoordinates)
{
	double norm_value=0;
	int n=numOfCoordinates;
	for(int i=0; i<totalPoints; i++)
	{
		for(int l=0; l<numOfCoordinates;l++)
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
