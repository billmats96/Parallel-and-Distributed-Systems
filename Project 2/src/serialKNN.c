/*
 * serialKNN.c
 *
 * Created on: Dec 11, 2017
 *      Author: linuxbill
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/time.h>
#include <math.h>


float* calculateSerialKnn(int k,float * array,int length,int numOfCoordinates);
float calculateDistance(int i, int j,float * array, int numOfCoordinates);
int* findNearestDistancesLabels(int k,float *array,float *nearestDistances, int length,int numOfCoordinates);
int cmpfunc (const void * a, const void * b);
void test(int k,int totalPoints,float *nearestDistances,int *labels);

struct timeval startwtime, endwtime;
double seq_time;
float *matlabDistances;
int *matlabLabels;


int main(int argc, char** argv)
{
	//Data will be read through a binary file named data.bin
	//Data is in float(32-bit) format
		FILE *fp;
		int sz,k,numOfCoordinates,length,totalPoints; //value of numOfCoordinates according to the coordinates of dataset to be used.
		float *dynarr, *nearestDistances;
		int *labels;
		int testMode=0; //set testMode 0 to use your own input.set testMode to 1 if you want to run tests using random inputs with matlab script.

		if(argc!=3){
			printf("Usage:./%s k n where k is the number of nearest neighbours and n the number of coordinates of every point in the dataset\n",argv[0]);
					  exit(1);
		}
		k=atoi(argv[1]);

		if(k>128)
		{
			printf("Too large value of k! Value should not exceed 128.");
			exit(1);
		}

		numOfCoordinates=atoi(argv[2]);

		//Opean data file.
		if(testMode)
		{
			fp=fopen("randInput.bin","rb");
		}
		else
		{

			//fp=fopen("/mnt/scratchdir/home/vmatsouk/totalFile.bin","rb");
			//fp=fopen("/mnt/scratchdir/home/vmatsouk/totalFile2.bin","rb");
			//fp=fopen("totalFile.bin","rb");
			fp=fopen("totalFile2.bin","rb");
		}

		if (!fp)
		{
			printf("Unable to open file!");
			return 1;
		}

		fseek(fp, 0, SEEK_END);
		sz = ftell(fp);
		printf("Size of data.bin in bytes=%d\n",sz);
		length= sz/sizeof(float); // take size in groups of floats
		printf("Length=%d\n",length);
		dynarr=(float *)malloc(length*sizeof(float)); //Dynamically allocate space in memory
		printf("Size of float = %d bytes\n",sizeof(float));

		if(dynarr==NULL)
		{
			printf("Cannot allocate memory..\n");
			exit(1);
		}

		rewind(fp);
		fread(dynarr,sizeof(float),length,fp);
		fclose(fp);

		//Print your input data

//		printf("\nThat's the data inserted:");
//		for (int i=0; i<(length) ; i++)
//		{
//			if(i%numOfCoordinates==0) printf("\n");
//			printf("%f ",dynarr[i]);
//		}
		printf("\n\n");
		totalPoints=length/numOfCoordinates;

		nearestDistances=(float *)malloc(k*totalPoints*sizeof(float));
		if(nearestDistances==NULL)
		{
			printf("Cannot allocate memory..\n");
			exit(1);
		}

		//Measure time of serial knn
		gettimeofday (&startwtime, NULL);
		nearestDistances=calculateSerialKnn(k,dynarr,length,numOfCoordinates);
		labels=(int *)(malloc(k*totalPoints*sizeof(int)));
		if(labels==NULL)
		{
			printf("Cannot allocate memory..\n");
			exit(1);
		}
		labels=findNearestDistancesLabels(k,dynarr,nearestDistances,length,numOfCoordinates);
		gettimeofday (&endwtime, NULL);
		seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6+ endwtime.tv_sec - startwtime.tv_sec);
		printf("Total time needed to execute serial knn: %f secs\n", seq_time);


		if(testMode)
		{
			matlabDistances=(float *)(malloc(k*totalPoints*sizeof(float)));
			matlabLabels=(int *)(malloc(k*totalPoints*sizeof(int)));

			fp=fopen("outputDISTtestKNN.bin","rb");
			fread(matlabDistances,sizeof(float),k*totalPoints,fp);
			fclose(fp);

			fp=fopen("outputIDXtestKNN.bin","rb");
			fread(matlabLabels,sizeof(int),k*totalPoints,fp);
			fclose(fp);
		}



		//Results. Uncomment to get them printed
//		printf("The distances i calculated:\n");
//		for(int i=0; i<(k*totalPoints); i++)
//		{
//			printf("Point %d is %f ",i/k,nearestDistances[i]); //prosthetoume +1 gia plhrh antistoixia me to matlab
//			printf("from point:%d\n",labels[i]);
//		}


		//Test serial knn.set testMode to true to run tests with matlab script
			 if(testMode)  test(k,totalPoints,nearestDistances,labels);


		//Free memory
		free((void *)labels);
		free((void *)nearestDistances);
		free((void *)dynarr);
		free((void *)matlabDistances);
		free((void *)matlabLabels);
		return 0;
}


int cmpfunc (const void * a, const void * b)
{
   return (  *(int*)a-*(int*)b );
}

//Calculates the distances of k nearest neighbours for each point

float* calculateSerialKnn(int k,float *array,int length,int numOfCoordinates)
{

	int totalPoints=length/numOfCoordinates;
	float temp=0;
	int counter=0;
	float *retDist; //holds the final distances from k nearest points, for each point
	float *tempRet; //holds the temporary distances of each point from all neighbours
	tempRet=(float *)(malloc((totalPoints*sizeof(float))));
	if(tempRet==NULL)
	{
		printf("Cannot allocate memory..\n");
		exit(1);
	}
	retDist=(float *)(malloc((k*totalPoints*sizeof(float))));
	if(retDist==NULL)
	{
		printf("Cannot allocate memory..\n");
		exit(1);
	}
	for(int i=0; i<totalPoints; i++)
	{
		counter=i*k;
		for(int j=0; j<totalPoints; j++)
		{
			if(i!=j)
			{
				temp=calculateDistance(i,j,array,numOfCoordinates);
				tempRet[j]=temp;
			}
			else
			{
				tempRet[j]=INFINITY; //1000000; Put a very large distance to signify distance from yourself
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

//Calculate distance between two points
float calculateDistance(int i, int j,float * array, int numOfCoordinates)
	{
		float distance=0;
		int n=numOfCoordinates;
		for(int l=0; l<numOfCoordinates;l++)
		{
			distance+=(array[n*i+l]-array[n*j+l])*(array[n*i+l]-array[n*j+l]);
		}
		//distance=sqrtf(distance);
		return distance;
	}
//Finds labels of k nearest points, for each point
int * findNearestDistancesLabels(int k,float *array,float *nearestDistances, int length,int numOfCoordinates)
{
	int totalPoints=length/numOfCoordinates;
	float temp=0;
	int counter=0,save;
	int *labelPointer;//holds the final labels from k nearest points, for each point
	labelPointer=(int *)(malloc(k*totalPoints*sizeof(int)));
	if(labelPointer==NULL)
	{
		printf("Cannot allocate memory..\n");
		exit(1);
	}
	for (int i=0; i<k*totalPoints; i++) labelPointer[i]=-1; //will help in sorting labels

	for (int i=0; i<totalPoints; i++)
	{
		counter=i*k; //scale according to the requested k
		for(int j=0; j<totalPoints; j++)
		{
			if(i!=j)
			{	//recalculate distance and search nearestDistances to identify that distance if it is eligible for neighbour.
				//save its label in the labelPointer array.
				temp=calculateDistance(i,j,array,numOfCoordinates);
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
				    	labelPointer[counter+save]=j;
				    }
			}
		}
	}
	return labelPointer;
}
void test(int k,int totalPoints,float *nearestDistances,int* labels)
{
   int flag=0;
	for(int i=0; i<k*totalPoints; i++)
   {
   	//printf("%d  and   %d\n",finLabels[i],matlabLabels[i]-1);
   	if(labels[i]!=matlabLabels[i]-1)
   	{
   		//printf("Spotted difference in labels. \n Matlab label: %d with distance:%1.6f\n "
   				//"Agorithm Label: %d with distance: %1.6f\n",matlabLabels[i]-1,matlabDistances[i],labels[i],nearestDistances[i]);

   		//printf("Deviation = %1.12f\n",fabsf(nearestDistances[i]-matlabDistances[i])); //have it print difference represented in 12 decimal digits.
   		    																//sometimes there's a small differnce to matlab and c floating point representation
   		    																//but the result is correct.
			if(fabsf(nearestDistances[i]-matlabDistances[i])>=0.00001)
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
      printf("Serial knn passed the test\n");
   else
	   printf("Serial knn failed the test\n");

}
