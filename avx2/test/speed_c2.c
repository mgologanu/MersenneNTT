#include<stdio.h>
#include<immintrin.h>
#include <stdalign.h>
#include "timing.h"
//#include <time.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

#include "ntt.h"

#define NR_RUNS 125*8*10000

//#define N 1048576

#define N 512
#define VAL_MAX 5000

int main()
{
  
  timing start;
  timing finish;
  timing t[2];

  
  double * result = (double *) malloc(NR_RUNS*sizeof(double));
  
  

  timing_now(&start);
  
   
  uint64_t *a = aligned_alloc(64, N*sizeof(uint64_t));
  
  int i=0;

  
  for (i=0; i<N; i++)
    {
      a[i] =  i;
    }
  

    
  for (i=0;  i<NR_RUNS; i++)
    {
      timing_now(&t[0]);
      
      v512(a);

      timing_now(&t[1]);
      result[i] = timing_diff(&t[1],&t[0]);
    }
  
  printf("%d ", N);
  for (i=0; i<3*8; i++)
    {
      printf(" %7.0f", result[i+NR_RUNS-3*8]);
    }
  printf("\n");

  timing_now(&finish);
  
  free(a);

    double mean = 0;
  double min = result[0];
  double max = result[0];

  int nb = 0;
  int nb2 = 0;
  
  for (i=0; i<NR_RUNS; i++)
    {
      if (result[i] < 0 || result[i] > VAL_MAX)
	{

	}
      else
	{
	  nb = nb + 1;
	  mean = mean +  result[i];
	  if (result[i] < min) min = result[i];
	  if (result[i] > max) max = result[i];
	}
    }
  mean = mean/nb;

  double std_sq = 0;
  for (i=0; i<NR_RUNS; i++)
    {
      if (result[i] < 0 || result[i] > VAL_MAX)
	{

	}
      else
	{
	  nb2 = nb2 + 1;
	  std_sq = std_sq + (result[i]-mean)*(result[i]-mean);
	}
    }
  
  printf("mean: %g, std: %g, min: %g, max: %g, nb: %d, nb2: %d\n",  mean, sqrt(std_sq/(nb2-1)), min, max, nb, nb2);


  for (i=0; i<NR_RUNS; i++)
    {
      if (result[i] < 0 || result[i] > VAL_MAX)
	{
	  
	}
      else
	{
	  //  printf("%g\n", result[i]);
	}
    }

  
}
