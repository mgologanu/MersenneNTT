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


#define NR_RUNS 125*8*10000

//#define NR_RUNS 1

#define N  NR_RUNS*8

#define PRIME_64 2147483647UL

#define PRIME_32 2147483647

#define VAL_MAX 100




static inline __attribute__((always_inline))
__m256i add_mod_mersenne_avx256_p31_32(__m256i a, __m256i b)
{

  __m256i t1, t2, t3, res;

  const __m256i p32 = _mm256_set1_epi32(PRIME_32);
  
  t3 = _mm256_add_epi32(a,b);
  
  t1 = _mm256_and_si256(t3, p32);

  t2 = _mm256_srli_epi32(t3, 31);
  
  t3  = _mm256_add_epi32(t1,t2);
    
  res = _mm256_and_si256(t3, p32);
  
  return res;
  
}


static inline __attribute__((always_inline))
__m256i sub_mod_mersenne_avx256_p31_32(__m256i a, __m256i b)
{

  __m256i t1, t2, t3, res;

  const __m256i p32 = _mm256_set1_epi32(PRIME_32);
    
  __m256i minus_b = _mm256_xor_si256(b, _mm256_set1_epi32(-1));
  
  t3 = _mm256_add_epi32(a,minus_b);
  
  t1 = _mm256_and_si256(t3, p32);

  t2 = _mm256_srli_epi32(t3, 31);
  
  t3  = _mm256_add_epi32(t1,t2);
    
  res = _mm256_and_si256(t3, p32);
  
  return res;
  
}



static inline __attribute__((always_inline))
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

  //const __m256i p64 = _mm256_set1_epi64x(PRIME_64);

  const __m256i p32 = _mm256_set1_epi32(PRIME_32);
  

  __m256i prod_even = _mm256_mul_epu32(a, b);
  
  // Move odd 32-bit lanes (1,3,5,7) down into the low half of each 64-bit lane
  __m256i a_odd = _mm256_srli_epi64(a, 32);
  __m256i b_odd = _mm256_srli_epi64(b, 32);
  
  __m256i prod_odd = _mm256_mul_epu32(a_odd, b_odd);
  
  __m256i hi_even = _mm256_srli_epi64(prod_even, 31);                       
  __m256i hi_odd  = _mm256_slli_epi64(_mm256_srli_epi64(prod_odd, 31), 32); 
    
  __m256i hi = _mm256_or_si256(hi_even, hi_odd);


  __m256i lo_odd = _mm256_slli_epi64(prod_odd, 32);
  
  const __m256i lo_mask = _mm256_setr_epi32(-1,  0, -1,  0, -1,  0, -1,  0); 

       
  __m256i lo1 = _mm256_or_si256(_mm256_and_si256(prod_even,lo_mask), lo_odd);

  __m256i lo =  _mm256_and_si256(lo1, p32);

  
  return  add_mod_mersenne_avx256_p31_32(lo, hi);
  

}


static inline __attribute__((always_inline))
void complex_mul_mod_mersenne_avx256_p31_32(__m256i a, __m256i b, __m256i c, __m256i s, __m256i * re, __m256i *im)
{

  const __m256i p32 = _mm256_set1_epi32(PRIME_32);
  
  __m256i ac = mul_mod_mersenne_avx256_p31_32(a, c);
  __m256i as = mul_mod_mersenne_avx256_p31_32(a, s);
  __m256i bc = mul_mod_mersenne_avx256_p31_32(b, c);
  __m256i bs = mul_mod_mersenne_avx256_p31_32(b, s);

  *re =  sub_mod_mersenne_avx256_p31_32(ac, bs);

  *im = add_mod_mersenne_avx256_p31_32(as, bc);

  return;
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
  uint32_t * d = aligned_alloc(32, N*sizeof(uint32_t));

  uint32_t * c = aligned_alloc(32, N*sizeof(uint32_t));
  uint32_t * s = aligned_alloc(32, N*sizeof(uint32_t));

  uint32_t * re = aligned_alloc(32, N*sizeof(uint32_t));
  uint32_t * im = aligned_alloc(32, N*sizeof(uint32_t));
 

  
  int32_t p=2147483647;


  for (int i= 0; i<N; i++)
    {
      a[i] = (p-1)-i;
      b[i] = i+1;
    }

  for (int i= 0; i<N; i++)
    {
      c[i] = (p-1)/2-i;
      s[i] = 2*i;
    }


  __m256i va, vb, vd, vc, vs, vre, vim;


  for (int i = 0; i < NR_RUNS; i++)
    //for (int i = 0; i < 1; i++)
  {

      

    //  t[i] = cpucycles();

    
        
    timing_now(&t[0]);

  
    va = _mm256_load_si256((__m256i *) &a[0+i*8]);
    
    vb = _mm256_load_si256((__m256i *) &b[0+i*8]);

    vc = _mm256_load_si256((__m256i *) &c[0+i*8]);

    vs = _mm256_load_si256((__m256i *) &s[0+i*8]);
    

    //vd = mul_mod_mersenne_avx256_p31_32(va, vb);
    //vc = mul_mod_mersenne_avx256_p31_32(va, vc);
    //vc = mul_mod_mersenne_avx256_p31_32(vb, vc);
    //vc = mul_mod_mersenne_avx256_p31_32(va, vc);

    complex_mul_mod_mersenne_avx256_p31_32(va, vb, vc, vs, &vre, &vim);
    





    
    _mm256_store_si256((__m256i *) &re[0+i*8], vre);
    _mm256_store_si256((__m256i *) &im[0+i*8], vim);

    timing_now(&t[1]);
    result[i] = timing_diff(&t[1],&t[0]);
    
  }

  /* printf("a = int64(["); */
  /* printf("%u", a[0]); */
  /* for (int i= 1; i<8; i++) */
  /*   { */
  /*     printf(", %u", a[i]); */
  /*   } */
  /* printf("]);\n"); */

  

  /* printf("b = int64(["); */
  /* printf("%u", b[0]); */
  /* for (int i= 1; i<8; i++) */
  /*   { */
  /*     printf(", %u", b[i]); */
  /*   } */
  /* printf("]);\n"); */
  


  /* printf("c = int64(["); */
  /* printf("%u", c[0]); */
  /* for (int i= 1; i<8; i++) */
  /*   { */
  /*     printf(", %u", c[i]); */
  /*   } */
  /* printf("]);\n"); */

  // print_results("Full convolution: ", t, N_TESTS);

  
     
 /*  for (int i= 0; i<8*N_TESTS; i++) */
 /*    { */
 /*        printf("%d, ", a[i]); */
 /*    } */
 /*  printf("\n"); */

 /*   for (int i= 0; i<8*N_TESTS; i++) */
 /*    { */
 /*        printf("%d, ", b[i]); */
 /*    } */
 /*  printf("\n"); */
  
 
  
 for (int i= 8*(NR_RUNS-1); i<8*NR_RUNS; i++)
    {
      printf("(%d, %d)  ", re[i], im[i]);
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


  for (int i=NR_RUNS-40; i<NR_RUNS; i++)
    {
      printf("%g ", result[i]);
    }

  printf("\n");


   /* for (int i=0; i<NR_RUNS; i++) */
  /*   { */
  /*     if (result[i] < 0 || result[i] > VAL_MAX) */
  /* 	{ */
	  
  /* 	} */
  /*     else */
  /* 	{ */
  /* 	  //  printf("%g\n", result[i]); */
  /* 	} */
  /*   } */

}

