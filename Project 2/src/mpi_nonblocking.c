/*
 * mpi_nonblocking.c
 *
 *  Created on: Dec 27, 2017
 *      Author: linuxbill
 */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>
#include "mpi.h"


#define MASTER 0


float* calculateSerialKnn(int k,float *array1,float *array2,int length,int numOfCoordinates,int offset2);
float calculateDistance(int i, int j,float * array1,float* array2, int numOfCoordinates);
int * findNearestDistancesLabels(int k,float *array1,float* array2,float *nearestDistances, int length,int numOfCoordinates,int offset,int offset2);
void update(float* tempDistances, int * tempLabels,int k,int totalPoints,int processes);
int cmpfunc (const void * a, const void * b);
void test(int rank,int k,int processes,int totalPoints);

struct timeval startwtime, endwtime;
double seq_time;
float *finDistances;
int *finLabels;
float *matlabDistances;
int *matlabLabels;

int main(int argc, char** argv)
{
	//Data will be read through a binary file named data.bin
	//Data is in float(32-bit) format

	//Variables declaration
	//value of numOfCoordinates according to the coordinates of dataset to be used.
	FILE *fp;
	float *dynarr,*dynarr_2,*receivebuf, *tempDistances;
	int numOfCoordinates,rank,size,value,offset,processes,length,totalPoints,chunk,offset2,sz,*tempLabels,k,len,prev,next,rec_value,tag1=1,tag2=2;
	char hostname[MPI_MAX_PROCESSOR_NAME];
	int store_value;
	double time,timecalc;
	MPI_Request reqs[4];   // required variable for non-blocking calls
	MPI_Status status[4];   // required variable for Waitall routine  
	
	int testMode=0; //set testMode 0 to use your own input.set testMode to 1 if you want to run tests using random inputs with matlab script.

	if (argc != 3) {
		printf("Usage:mpiexec -np p %s k n  where k is the number of nearest neighbours and n the number of coordinates of every point in the dataset\n",argv[0]);
		  exit(1);
	  }

	k=atoi(argv[1]);
	numOfCoordinates=atoi(argv[2]);
	//value=processes-1;

	if(k>128)
	{
		printf("Too large value of k! Value should not exceed 128.");
		exit(1);
	}
	//Initialize Message Passing Interface environment
	MPI_Init (&argc,&argv);
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Get_processor_name(hostname, &len);

	processes=size; //just for convenience
	value=size-1;

	//Opean data file.
	if(testMode)
	{
		fp=fopen("randInput.bin","rb");
	}
	else
	{
		//fp=fopen("/mnt/scratchdir/home/vmatsouk/totalFile.bin","rb");
		//fp=fopen("/mnt/scratchdir/home/vmatsouk/totalFile2.bin","rb");
		fp=fopen("totalFile.bin","rb");
		//fp=fopen("totalFile2.bin","rb");
	}

	if (!fp)
	{
		printf("Unable to open file!");
		return 1;
	}
	fseek(fp, 0, SEEK_END);
	sz = ftell(fp); //obtain file size in bytes
	//printf("Size of data.bin in bytes=%d\n",sz);
	length= sz/sizeof(float); // take size in groups of floats
	//printf("Length=%d\n",length);
	totalPoints=length/numOfCoordinates;
	if(totalPoints%processes!=0)
	{
		printf("Error! Number of points must be divisible by processes.");
		exit(1);
	}

	printf ("Number of tasks= %d My rank= %d Running on %s\n", size,rank,hostname);
	chunk=length/processes;//the size of block that will be distributed to the process
	dynarr=(float *)malloc(chunk*sizeof(float)); //Data to be held by each rank
	dynarr_2=(float *)malloc(chunk*sizeof(float)); //Memory allocation for the distributed data
	receivebuf=(float *)malloc(chunk*sizeof(float));//Buffer for the incoming dataset
	//printf("Size of float = %d bytes\n",sizeof(float));

	if(dynarr==NULL)
	{
		printf("Cannot allocate memory..\n");
		exit(1);
	}
	rewind(fp);
	fseek(fp,rank*chunk*sizeof(float),SEEK_SET);  //Each rank reads and stores the corresponding part of file
	fread(dynarr,sizeof(float),chunk,fp);
	rewind(fp);
	fseek(fp,rank*chunk*sizeof(float),SEEK_SET);  //Initially the distributed data is the same block to execute knn with itself.
	fread(dynarr_2,sizeof(float),chunk,fp);
	fclose(fp);

	//Input data for each process
	//	printf("\nThat's the data block held by process %d:", rank);
	//	for (int i=0; i<(chunk) ; i++)
	//	{
	//		if(i%numOfCoordinates==0) printf("\n");
	//		printf("%f ",dynarr[i]);
	//	}
	//printf("\n\n");

	tempLabels=(int *)(malloc(k*totalPoints/processes*sizeof(int)));//final distances array for each process
	tempDistances=(float *)malloc(k*totalPoints/processes*sizeof(float));//final labels array for each process
	finDistances=(float *)malloc(k*totalPoints/processes*sizeof(float));//temporary distances array for each process
	finLabels=(int *)(malloc(k*totalPoints/processes*sizeof(int)));//temporary labels array for each process

	if(testMode)
	{
		matlabDistances=(float *)(malloc(k*totalPoints/processes*sizeof(float)));
		matlabLabels=(int *)(malloc(k*totalPoints/processes*sizeof(int)));


		fp=fopen("outputDISTtestKNN.bin","rb");
		fseek(fp,rank*k*totalPoints/processes*sizeof(float),SEEK_SET);
		fread(matlabDistances,sizeof(float),k*totalPoints/processes,fp);
		fclose(fp);

		fp=fopen("outputIDXtestKNN.bin","rb");
		fseek(fp,rank*k*totalPoints/processes*sizeof(int),SEEK_SET);
		fread(matlabLabels,sizeof(int),k*totalPoints/processes,fp);
		fclose(fp);
	}


	int blocks[processes]; //a simple array that will help in offset and offset2 calculation
	for (int i=0; i<processes;i++)
			blocks[i]=i*chunk;

	for(int j=0; j<k*totalPoints/processes; j++) finDistances[j]=1000000; //Initialization that will help in sorting. final distances in our
																			//application will always be much lower.
	
	//Synchronise processes and start stopwatch.
	MPI_Barrier (MPI_COMM_WORLD);
	if (rank==MASTER){
		gettimeofday (&startwtime, NULL);
		timecalc-=MPI_Wtime();
	}

	//Determine left and right neighbors
	prev = rank-1;
	next = rank+1;
	if (rank == MASTER)  prev = size - 1;
	if (rank == (size - 1))  next = 0;

	//Rec_value is the block number which is being passed around by the processes. Initialized here.
	rec_value=(prev+1)%size;

	 MPI_Isend(&rec_value, 1, MPI_INT, next, tag1, MPI_COMM_WORLD, &reqs[2]);
	 MPI_Isend(dynarr_2,chunk,MPI_FLOAT,next,tag2,MPI_COMM_WORLD, &reqs[3]);
	 MPI_Irecv(&store_value, 1, MPI_INT, prev, tag1, MPI_COMM_WORLD, &reqs[0]);
	 MPI_Irecv(receivebuf,chunk,MPI_FLOAT,prev,tag2,MPI_COMM_WORLD,&reqs[1]);

	 //execute knn for each process with itself
	 //offset variables help us specify the blocks of data we're working with
	offset=blocks[rank]/numOfCoordinates;
	offset2=0;
	finDistances=calculateSerialKnn(k,dynarr,dynarr_2,chunk,numOfCoordinates,offset2);
	finLabels=findNearestDistancesLabels(k,dynarr,dynarr_2,finDistances,chunk,numOfCoordinates,offset,offset2);
	//update(tempDistances,tempLabels,k,totalPoints,processes);
	MPI_Waitall(4, reqs, status);

	for(int i=0; i<chunk; i++)
	{
		dynarr_2[i]=receivebuf[i];
	}
	rec_value=store_value;

	//Wait again for all processes to synchronise
	 MPI_Barrier (MPI_COMM_WORLD);
	 if(rank==0) timecalc+=MPI_Wtime(); //take into account the initial calculations.

	 while(value>0){
		 MPI_Barrier (MPI_COMM_WORLD);
	   MPI_Isend(&rec_value, 1, MPI_INT, next, tag1, MPI_COMM_WORLD, &reqs[2]);
	   MPI_Isend(dynarr_2,chunk,MPI_FLOAT,next,tag2,MPI_COMM_WORLD, &reqs[3]);
	   MPI_Irecv(&store_value, 1, MPI_INT, prev, tag1, MPI_COMM_WORLD, &reqs[0]);
	   MPI_Irecv(receivebuf,chunk,MPI_FLOAT,prev,tag2,MPI_COMM_WORLD,&reqs[1]);

	   MPI_Barrier (MPI_COMM_WORLD);
	   if(rank==0) timecalc-=MPI_Wtime();

	   offset=blocks[rec_value]/numOfCoordinates;
	   offset2=blocks[abs(rec_value-rank)]/numOfCoordinates;
	   tempDistances=calculateSerialKnn(k,dynarr,dynarr_2,chunk,numOfCoordinates,offset2);
	   tempLabels=findNearestDistancesLabels(k,dynarr,dynarr_2,tempDistances,chunk,numOfCoordinates,offset,offset2);
	   update(tempDistances,tempLabels,k,totalPoints,processes);
	   value--;
	   MPI_Waitall(4, reqs, status);

	  for(int i=0; i<chunk; i++)
	  {
		 dynarr_2[i]=receivebuf[i];
	  }
		 rec_value=store_value;

	  MPI_Barrier (MPI_COMM_WORLD);
	  if(rank==0) timecalc+=MPI_Wtime();
	 }
	 
 //Continue - do more work
	MPI_Barrier (MPI_COMM_WORLD);
	if (rank==MASTER)
	{
		gettimeofday (&endwtime, NULL);
		seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6+ endwtime.tv_sec - startwtime.tv_sec);
			printf("Total time needed to execute parallel knn with mpi and %d processes: %f secs\n",size, seq_time);
			printf("Calculations time=%f\n",timecalc);
	}
	MPI_Finalize();

  // Each process prints its results.
//	printf("\nThe distances calculated by process %d:\n",rank);
//	int counter=0;
//	for(int j=0; j<k*totalPoints/processes; j++)
//	{
//		printf("\nMin distance for point %d= %f from point %d",counter/k+rank*totalPoints/processes+1,finDistances[j],finLabels[j]+1);
//		counter++;
//	}
//	printf("\n");

//Each process tests its results.set testMode to true to run tests with matlab script
  if(testMode)  test(rank,k,processes,totalPoints);

  //Each process frees allocated memory

  	free((void *)tempDistances);
  	free((void *)tempLabels);
  	free((void *)finDistances);
  	free((void *)finLabels);
  	free((void *)matlabDistances);
  	free((void *)matlabLabels);
  	free((void *)dynarr);
  	free((void *)dynarr_2);
  	free((void *)receivebuf);
  	return 0;
 }


 int cmpfunc (const void * a, const void * b)
 {
	return (  *(int*)a-*(int*)b );
 }
 //For each point in array1, calculates serially the distances of k nearest neighbours from array2
 float* calculateSerialKnn(int k,float *array1,float *array2,int length,int numOfCoordinates,int offset2)
 {
	int totalPoints=length/numOfCoordinates;
	float temp=0;
	int counter=0;
	float *retDist;
	float *tempRet;
	tempRet=(float *)(malloc((totalPoints*sizeof(float))));//holds the final distances from k nearest points, for each point
	retDist=(float *)(malloc((k*totalPoints*sizeof(float))));//holds the temporary distances of each point from all neighbours

	for(int i=0; i<totalPoints; i++)
	{
		counter=i*k;

		for(int j=0; j<totalPoints; j++)
		{
			if(i!=(j+offset2))
			{
				temp=calculateDistance(i,j,array1,array2,numOfCoordinates);
				tempRet[j]=temp;
			}
			else
			{
				tempRet[j]=INFINITY; // 1000000; Put a very large distance to signify distance from yourself
									//and prevent from being identied as neighbour.Could be replaced by value=INFINITY
			}
		}
		//Sort calculated distances from all its neighbours and then save the k-nearest
		qsort((void *)&tempRet[0],totalPoints,sizeof(float),cmpfunc);
		for(int j=0; j<k;j++)
		{
			retDist[counter+j]=tempRet[j];
		}
	}
	free(tempRet);
	return retDist;
 }

 //Calculate distance between two points from array1 array2
 float calculateDistance(int i, int j,float * array1,float* array2, int numOfCoordinates)
	{
		float distance=0;
		int n=numOfCoordinates;
		for(int l=0; l<numOfCoordinates;l++)
		{
			distance+=(array1[n*i+l]-array2[n*j+l])*(array1[n*i+l]-array2[n*j+l]);
		}
		return distance;
	}

 //For each point in array1,finds labels of k nearest points from array2
 int * findNearestDistancesLabels(int k,float *array1,float* array2,float *nearestDistances, int length,int numOfCoordinates,int offset,int offset2)
 {
	int totalPoints=length/numOfCoordinates;
	float temp=0;
	int counter=0,save;
	int *labelPointer;
	labelPointer=(int *)(malloc(k*totalPoints*sizeof(int)));//holds the final labels from k nearest points, for each point
	for (int i=0; i<k*totalPoints; i++) labelPointer[i]=-1; //will help in sorting labels

	for (int i=0; i<totalPoints; i++)
	{
		counter=i*k;//scale according to the requested k

		for(int j=0; j<totalPoints; j++)
		{
			if(i!=(j+offset2))
			{
				temp=calculateDistance(i,j,array1,array2,numOfCoordinates);//recalculate distance and search nearestDistances to identify that distance if it is eligible for neighbour.
																			//save its label in the labelPointer array.
				if (temp<=nearestDistances[counter+k-1])
				 {
					for(int z=0; z<=(k-1); z++)
					{
						 if(temp==nearestDistances[counter+z])
						 {
							save=z;
							break;
						 }
					}
					while (labelPointer[counter+save]>0 && save<=k-2) //<=k-2 cause we want save++ to go till k-1
						save++;
				labelPointer[counter+save]=j+offset;
				}
			}
		}
	}
	return labelPointer;
 }

 //Compares tempDistances and tempLabels with the so-far k nearest neighbours data and updates finDistances and finLabels arrays
 void update(float* tempDistances, int* tempLabels,int k,int totalPoints,int processes)
 {
	float min;
	int label;
	int counter=0,indx=0,extra_counter=-1;

	min=finDistances[k-1];
	label=finLabels[k-1];

	for(int j=0; j<k*totalPoints/processes; j++)
	{
		extra_counter++;
		if(extra_counter==k)
		{
			min=finDistances[j+k-1];
			label=finLabels[j+k-1];
			counter++;
			extra_counter=0;
		}
		if (tempDistances[j]<=min) //if we use <=, we take the furthest label. if we use <,we take the first label
		{
			for(int i=(counter*k); i<=((counter+1)*k-1); i++)
			{
				if(tempDistances[j]<=finDistances[i])
				{
					indx=i;
					break;
				}
			}

			for(int i=((counter+1)*k-1); i>indx; i--)
			{
				finLabels[i]=finLabels[i-1];
				finDistances[i]=finDistances[i-1];
			}
			finDistances[indx]=tempDistances[j];
			finLabels[indx]=tempLabels[j];

			min= finDistances[(counter+1)*k-1];
			label=finLabels[(counter+1)*k-1];

		}
	}
 }

 void test(int rank,int k,int processes,int totalPoints)
 {
	   int flag=0;
		for(int i=0; i<k*totalPoints/processes; i++)
	    {
	    	//printf("%d  and   %d\n",finLabels[i],matlabLabels[i]-1);
	    	if(finLabels[i]!=matlabLabels[i]-1)
	    	{
	    		//printf("Spotted difference in labels. \n Matlab label: %d with distance:%1.6f\n "
	    				//"Agorithm Label: %d with distance: %1.6f\n",matlabLabels[i]-1,matlabDistances[i],finLabels[i],finDistances[i]);

	    		//printf("Deviation = %1.12f\n",fabsf(finDistances[i]-matlabDistances[i])); //have it print difference represented in 12 decimal digits.
	    	   																//sometimes there's a small differnce to matlab and c floating point representation
	    		    																//but the result is correct.
				if(fabsf(finDistances[i]-matlabDistances[i])>=0.000001)
				{
					printf("Considerable error in distances.\n");
					flag=1;
					break;
				}
				else
				{
					//printf("it's ok\n");
				}
	    	}

	    }

	    if(!flag)
	      printf("Process %d passed the test\n",rank);
	   else
		   printf("Process %d failed the test\n",rank);


 }
	
