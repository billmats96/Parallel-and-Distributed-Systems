/*
 * pthreadbitonic_imp.c
 *
 *Author: Matsoukas Vasileios
 *        Undergraduate Student, Department of Electrical and Computer Engineering
 *		   Aristotle University of Thessaloniki, Greece
 * AEM:8743
 * email: vmatsouk@auth.gr
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>
struct timeval startwtime, endwtime;
double seq_time;


int N;          // data array size
int *a;         // data array to be sorted
int NTHREADS;   // threads to be created
const int ASCENDING  = 1;
const int DESCENDING = 0;


void init(void);
void print(void);
void test(void);
inline void exchange(int i, int j);
void compare(int i, int j, int dir);
void impBitonicSort(int);

//a struct that contains data about a thread
struct thread_data{
   int  tid; 				//thread id
   int  iterations ;		//iterations to execute by thread
   int  k;					//k variable used by thread in for loop
   int  j;					//j variable used by thread in for loop
};


void *do_work(void *threadarg)
{
  int i, start, mytid, end,j,k,ITER;
  struct thread_data *my_data;
  my_data = (struct thread_data *) threadarg; //get a pointer to the current thread
  mytid=my_data->tid;   					  //obtain thread's ID
  ITER=my_data->iterations;					  //iterations to be done by current thread
  start = mytid *ITER;						  //start iterations from...
  end = start + ITER;						  //end of iterations
  j=my_data->j;
  k=my_data->k;
  //printf ("Thread %d doing iterations %d to %d\n",mytid,start,end-1);
  for (i=start; i<end; i++) {
  	int ij=i^j;
  	if ((ij)>i) {
  	  if ((i&k)==0 && a[i] > a[ij])
  	      exchange(i,ij);
  	  if ((i&k)!=0 && a[i] < a[ij])
  	      exchange(i,ij);
  	           }
        }
  pthread_exit(NULL);
}


/** the main program **/
int main(int argc, char **argv) {

  if (argc != 3) {
	printf("Usage: %s q p\n  where n=2^q is problem size (power of two) and t=2^p the number of threads\n",
			   argv[0]);
    exit(1);
  }
  N = 1<<atoi(argv[1]);
  NTHREADS=1<<atoi(argv[2]);
  a = (int *) malloc(N * sizeof(int));
  printf("Threads:%d\n",NTHREADS);

  init();
  gettimeofday (&startwtime, NULL);
  impBitonicSort(NTHREADS);
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

 printf("Imperative wall clock time = %f\n", seq_time);

  test();
 // printf("%f", seq_time);

  // print();
}

/** -------------- SUB-PROCEDURES  ----------------- **/

/** procedure test() : verify sort results **/
void test() {
  int pass = 1;
  int i;
  for (i = 1; i < N; i++) {
    pass &= (a[i-1] <= a[i]);
  }

  printf(" TEST %s\n",(pass) ? "PASSed" : "FAILed");
}


/** procedure init() : initialize array "a" with data **/
void init() {
  int i;
  for (i = 0; i < N; i++) {
    a[i] = rand() % N; // (N - i);
  }
}

/** procedure  print() : print array elements **/
void print() {
  int i;
  for (i = 0; i < N; i++) {
    printf("%d\n", a[i]);
  }
  printf("\n");
}


/** INLINE procedure exchange() : pair swap **/
inline void exchange(int i, int j) {
  int t;
  t = a[i];
  a[i] = a[j];
  a[j] = t;
}



/** procedure compare()
   The parameter dir indicates the sorting direction, ASCENDING
   or DESCENDING; if (a[i] > a[j]) agrees with the direction,
   then a[i] and a[j] are interchanged.
**/
inline void compare(int i, int j, int dir) {
  if (dir==(a[i]>a[j]))
    exchange(i,j);
}




/*
  imperative version of bitonic sort
*/
void impBitonicSort(int NTHREADS) {

  int j,k;
  int ITERATIONS=N/NTHREADS; //iterations will be shared equally by threads
  struct thread_data thread_data_array[NTHREADS];
  for (k=2; k<=N; k=2*k) {
    for (j=k>>1; j>0; j=j>>1) {
    	  int i, start;
    	  pthread_t threads[NTHREADS];
    	  pthread_attr_t attr;

    	  pthread_attr_init(&attr);
    	  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE); //set thread's attribute to joinable
    	  for(i=0; i<NTHREADS; i++)
    	  {
    		  thread_data_array[i].tid=i;
    		  thread_data_array[i].j = j;
    		  thread_data_array[i].k = k;
    		  thread_data_array[i].iterations=ITERATIONS;
    	  }
    	  for (i=0; i<NTHREADS; i++)
    	  {
    	    pthread_create(&threads[i], &attr, do_work, (void *) &thread_data_array[i]);
    	  }
    	  for (i=0; i<NTHREADS; i++)
    	  {
    		 // printf("I am thread: %d\n",thread_data_array[i].tid);
    	    pthread_join(threads[i], NULL);
    	  }
    }
  }
}


