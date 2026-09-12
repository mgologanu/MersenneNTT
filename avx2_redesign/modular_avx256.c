#include<stdint.h>
#include <stdio.h>


#include<immintrin.h>
#include <stdalign.h>



#include "cpucycles.h"
#include "speed_print.h"

#define N_TESTS 100000000

uint64_t t[N_TESTS];



#define PRIME_64 2147483647UL

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

  //_Alignas(32) int32_t a[8*N_TESTS];
  //_Alignas(32) int32_t b[8*N_TESTS];
  //_Alignas(32) int32_t c[8*N_TESTS];

  uint32_t * a = aligned_alloc(32, 8*N_TESTS*sizeof(uint32_t));
  uint32_t * b = aligned_alloc(32, 8*N_TESTS*sizeof(uint32_t));
  uint32_t * c = aligned_alloc(32, 8*N_TESTS*sizeof(uint32_t));

  int32_t *pa, *pb, *pc;
  
  int64_t cc[4];

  int32_t p=2147483647;


  for (int i= 0; i<8*N_TESTS; i++)
    {
      a[i] = (p-1)-i;
      b[i] = i+1;
    }


  __m256i va, vb, vc;


  pa = a;
  pb = b;
  pc = c;

  for (int i = 0; i < N_TESTS; i++)
  {

    //Start measuring speed
      

    //  t[i] = cpucycles();

  
      
    va = _mm256_load_si256((__m256i *) &a[0+i*8]);
    
    vb = _mm256_load_si256((__m256i *) &b[0+i*8]);
    
    
    vc = mul_mod_mersenne_avx256_p31_32(va, vb);
    
    _mm256_store_si256((__m256i *) &c[0+i*8], vc);
    
    //    pa = pa + 8;
    // pb = pb + 8;
    // pc = pc + 8;
    
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
  
  
 for (int i= 8*(N_TESTS-1); i<8*N_TESTS; i++)
    {
        printf("%d, ", c[i]);
    }
  printf("\n");

  
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
