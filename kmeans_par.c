//December 2020
//Marcelo dos Santos
//marcelouepg@gmail.com

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#define DIM 3
#define NUM_THREADS 4//number of threads used on kmeans

#include <stdlib.h>
#include <sys/time.h>

int main(void) {

	int i, j, k, n;
	double *x, *mean, *sum;
	int *cluster, *count, *flip_par;
	int flips,ret;
	int num_threads;
	ret=scanf("%d", &k);
	ret=scanf("%d", &n);
	x = (double *)malloc(sizeof(double)*DIM*n);
	mean = (double *)malloc(sizeof(double)*DIM*k);
	sum= (double *)malloc(sizeof(double)*DIM*k);
	cluster = (int *)malloc(sizeof(int)*n);
	count = (int *)malloc(sizeof(int)*k);
	flip_par = (int *)malloc(sizeof(int)*NUM_THREADS);
	for (i = 0; i<k; i++)
		ret=scanf("%lf %lf %lf", mean+i*DIM, mean+i*DIM+1, mean+i*DIM+2);
	for (i = 0; i<n; i++){
		ret=scanf("%lf %lf %lf", x+i*DIM, x+i*DIM+1, x+i*DIM+2);
		cluster[i] = 0;
	}
	for (i=0;i<NUM_THREADS;i++){flip_par[i]=0;}
	flips = n;
	
	while (flips>0) {

		flips = 0;
		for (j = 0; j < k; j++) {
			count[j] = 0;
			for (i = 0; i < DIM; i++)
				sum[j*DIM+i] = 0.0;
		}

		//parallel region starts here
		#pragma omp parallel num_threads(NUM_THREADS)
		{
			int nthreads=omp_get_num_threads();
			int ID=omp_get_thread_num();
			int tam_p_thread=n/nthreads;
			int min, max,ip,colorp,cp,jp;
			double dmin, dx;
			min=tam_p_thread*ID;
			if(ID!=nthreads-1){max=tam_p_thread*(ID+1);}
			else{max=n;}
			if(ID==0){num_threads=nthreads;}
			
			for (ip = min; ip < max; ip++) {
				dmin = -1; colorp = cluster[ip];
				for (cp = 0; cp < k; cp++) {
					dx = 0.0;
					for (jp = 0; jp < DIM; jp++)
						dx +=  (x[ip*DIM+jp] - mean[cp*DIM+jp])*(x[ip*DIM+jp] - mean[cp*DIM+jp]);
					if (dx < dmin || dmin == -1) {
						colorp = cp;
						dmin = dx;
					}
				}
				if (cluster[ip] != colorp) {
					flip_par[ID]++;
					cluster[ip] = colorp;
			  		}
			}
			
		}
		//parallel region ends here
		for (i=0;i<num_threads;i++){flips=flips+flip_par[i];flip_par[i]=0;}

		for (i = 0; i < n; i++) {
			count[cluster[i]]++;
			for (j = 0; j < DIM; j++){
				sum[cluster[i]*DIM+j] += x[i*DIM+j];

			}
		}

		for (i = 0; i < k; i++) {
			for (j = 0; j < DIM; j++) {
				mean[i*DIM+j] = sum[i*DIM+j]/count[i];
  			}
		}

	}
	ret=ret+0;

	for (i = 0; i < k; i++) {
		for (j = 0; j < DIM; j++)
			printf("%5.2f ", mean[i*DIM+j]);
		printf("\n");
	}
	#ifdef DEBUG
	for (i = 0; i < n; i++) {
		for (j = 0; j < DIM; j++)
			printf("%5.2f ", x[i*DIM+j]);
		printf("%d\n", cluster[i]);
	}
	#endif

	return(0);
}
