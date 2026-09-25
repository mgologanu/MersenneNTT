#ifndef _MRSN_AVX2_INTERNAL_H
#define _MRSN_AVX2_INTERNAL_H


#define MRSN_PRIME 2147483647

#define MRSN_PRIME_POW 31



static inline __attribute__((always_inline))
__m256i add(__m256i a, __m256i b)
{

  __m256i t1, t2, t3, res;

  const __m256i p = _mm256_set1_epi32(MRSN_PRIME);
  
  t3 = _mm256_add_epi32(a, b);
  
  t1 = _mm256_and_si256(t3, p);

  t2 = _mm256_srli_epi32(t3, MRSN_PRIME_POW);
  
  t3 = _mm256_add_epi32(t1, t2);
    
  res = _mm256_and_si256(t3, p);
  
  return res;
  
}



static inline __attribute__((always_inline))
__m256i sub(__m256i a, __m256i b)
{

  __m256i t1, t2, t3, res;

  const __m256i p = _mm256_set1_epi32(MRSN_PRIME);
    
  __m256i minus_b = _mm256_xor_si256(b, _mm256_set1_epi32(MRSN_PRIME));
  
  t3 = _mm256_add_epi32(a, minus_b);
  
  t1 = _mm256_and_si256(t3, p);

  t2 = _mm256_srli_epi32(t3, MRSN_PRIME_POW);
  
  t3  = _mm256_add_epi32(t1, t2);
    
  res = _mm256_and_si256(t3, p);
  
  return res;
  
}




static inline __attribute__((always_inline))
__m256i mul(__m256i a, __m256i b)
{
  const __m256i p = _mm256_set1_epi32(MRSN_PRIME);
  
  const __m256i lo_mask = _mm256_setr_epi32(-1,  0, -1,  0, -1,  0, -1,  0); 

    
  __m256i prod_even = _mm256_mul_epu32(a, b);
  
  __m256i a_odd = _mm256_srli_epi64(a, 32);
  __m256i b_odd = _mm256_srli_epi64(b, 32);
  
  __m256i prod_odd = _mm256_mul_epu32(a_odd, b_odd);
  
  __m256i hi_even = _mm256_srli_epi64(prod_even, MRSN_PRIME_POW);                       
  __m256i hi_odd  = _mm256_slli_epi64(_mm256_srli_epi64(prod_odd, MRSN_PRIME_POW), 32); 
  
  __m256i hi = _mm256_or_si256(hi_even, hi_odd);
  
  __m256i lo_odd = _mm256_slli_epi64(prod_odd, 32);
  
  
  __m256i lo1 = _mm256_or_si256(_mm256_and_si256(prod_even,lo_mask), lo_odd);
  
  __m256i lo =  _mm256_and_si256(lo1, p);
  
  
  return  add(lo, hi);
  
}




static inline __attribute__((always_inline))
__m256i rot15(__m256i a)
{
  //   return ((x << 15) | (x >> 16)) & 0x7FFFFFFF;

  __m256i t1, t2, t3, res;

  const __m256i p = _mm256_set1_epi32(MRSN_PRIME);
   
  t1 = _mm256_slli_epi32(a, 15);

  t2 = _mm256_srli_epi32(a, 16);

  t3 = _mm256_or_si256(t1, t2);

  res =  _mm256_and_si256(t3, p);

  return res;
}



static inline __attribute__((always_inline))
__m256i rot(__m256i a, int k)
{
  //   return ((x << k) | (x >> 31-k)) & 0x7FFFFFFF;

  __m256i t1, t2, t3, res;

  int km =  MRSN_PRIME_POW - k;

  const __m256i p = _mm256_set1_epi32(MRSN_PRIME);
   
  t1 = _mm256_slli_epi32(a, k);

  t2 = _mm256_srli_epi32(a, km);

  t3 = _mm256_or_si256(t1, t2);

  res =  _mm256_and_si256(t3, p);

  return res;
}

#endif


static inline __attribute__((always_inline))
bf1(__m256i * a_re, __m256i * a_im, __m256i * b_re, __m256i * b_im)
{
  /*
    a = a + b
    b = a - b
  */

  __m256i t1, t2;
  
  t1    = add(*a_re, *b_re);
  *b_re = sub(*a_re, *b_re);
  *a_re = t1;

  t2    = add(*a_im, *b_im);
  *b_im = sub(*a_im, *b_im);
  *a_im = t2;

}

