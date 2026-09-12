#include<stdint.h>
#include <stdio.h>
#include <math.h>

#include<immintrin.h>
#include <stdalign.h>

#include "timing.h"


/* #include "cpucycles.h" */
/* #include "speed_print.h" */

/* #define N_TESTS 100000000 */

/* uint64_t t[N_TESTS]; */


#define NR_RUNS 125*8*100000

#define N  NR_RUNS*8

#define PRIME_64 2147483647UL

#define VAL_MAX 100

__m256i mul_mod_mersenne_avx256_p31_32(__m256i a, __m256i b)
{
  /*

    Input: a and b unsigned reminders mod 2^31-1

    0 <= a b <= 2^31-1
    b <= a b <= 2^31-1

    Note: Accepting double zero where a and/or b = p = 2^31-1.
    
    Output: Unsigned reminder mod 2^31-1


    Widening multiplication, followed by 2 additions for a full reduction
    
  */    

  const __m256i p64 = _mm256_set1_epi64x(PRIME_64);
  

  __m256i prod_even = _mm256_mul_epu32(a, b);
  
  // Move odd 32-bit lanes (1,3,5,7) down into the low half of each 64-bit lane
  __m256i a_odd = _mm256_srli_epi64(a, 32);
  __m256i b_odd = _mm256_srli_epi64(b, 32);
  
  __m256i prod_odd = _mm256_mul_epu32(a_odd, b_odd);
  
  
  
  __m256i t1_even = _mm256_and_si256(prod_even, p64);
  __m256i t2_even = _mm256_srli_epi64(prod_even, 31);
  __m256i t3_even = _mm256_add_epi64(t1_even,t2_even);

   t1_even = _mm256_and_si256(t3_even, p64);
   t2_even = _mm256_srli_epi64(t3_even, 31);

   t3_even = _mm256_add_epi64(t1_even,t2_even);

     
  __m256i t1_odd = _mm256_and_si256(prod_odd, p64);
  __m256i t2_odd = _mm256_srli_epi64(prod_odd, 31);
  __m256i t3_odd = _mm256_add_epi64(t1_odd,t2_odd);

   t1_odd = _mm256_and_si256(t3_odd, p64);
   t2_odd = _mm256_srli_epi64(t3_odd, 31);

   t3_odd = _mm256_add_epi64(t1_odd,t2_odd);

   __m256i prod_mod_even =  _mm256_and_si256(t3_even,p64);
   
   __m256i prod_mod_odd  = _mm256_slli_epi64(_mm256_and_si256(t3_odd, p64), 32); // lanes 1,3,5,7

   __m256i  res = _mm256_or_si256(prod_mod_even, prod_mod_odd);
    

    return res;
  

}


int main()
{


  
  timing start;
  timing finish;
  timing t[2];

  
  double * result = (double *) malloc(NR_RUNS*sizeof(double));
  
  

  
  //_Alignas(32) int32_t a[8*N_TESTS];
  //_Alignas(32) int32_t b[8*N_TESTS];
  //_Alignas(32) int32_t c[8*N_TESTS];

  uint32_t * a = aligned_alloc(32, N*sizeof(uint32_t));
  uint32_t * b = aligned_alloc(32, N*sizeof(uint32_t));
  uint32_t * c = aligned_alloc(32, N*sizeof(uint32_t));

  int32_t p=2147483647;


  for (int i= 0; i<N; i++)
    {
      a[i] = (p-1)-i;
      b[i] = i+1;
    }


  __m256i va, vb, vc;


  for (int i = 0; i < NR_RUNS; i++)
  {

      

    //  t[i] = cpucycles();

    va = _mm256_load_si256((__m256i *) &a[0+i*8]);
    
    vb = _mm256_load_si256((__m256i *) &b[0+i*8]);
    
    timing_now(&t[0]);
    vc = mul_mod_mersenne_avx256_p31_32(va, vb);
    timing_now(&t[1]);
    result[i] = timing_diff(&t[1],&t[0]);
    
    _mm256_store_si256((__m256i *) &c[0+i*8], vc);

   
    
  }

 




  // print_results("Full convolution: ", t, N_TESTS);

  
     
  /* for (int i= 0; i<8*N_TESTS; i++) */
  /*   { */
  /*       printf("%d, ", a[i]); */
  /*   } */
  /* printf("\n"); */

  /*  for (int i= 0; i<8*N_TESTS; i++) */
  /*   { */
  /*       printf("%d, ", b[i]); */
  /*   } */
  /* printf("\n"); */
  
  
 for (int i= 8*(NR_RUNS-1); i<8*NR_RUNS; i++)
    {
        printf("%d, ", c[i]);
    }
  printf("\n");

  
  double mean = 0;
  double min = result[0];
  double max = result[0];

  int nb = 0;
  int nb2 = 0;
  
  for (int i=0; i<NR_RUNS; i++)
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
  for (int i=0; i<NR_RUNS; i++)
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


  for (int i=0; i<NR_RUNS; i++)
    {
      if (result[i] < 0 || result[i] > VAL_MAX)
	{
	  
	}
      else
	{
	  //  printf("%g\n", result[i]);
	}
    }


  
  /* for (int i= 0; i<N_TESTS; i++) */
  /*   { */
  /*     printf("%lu, ", t[i]); */
  /*   } */
  /* printf("\n"); */
  
  
  /*
  for (int i = 0; i<4; i++)
    printf("%ld, ", cc[i]);

  printf("\n");
  */
  
}
