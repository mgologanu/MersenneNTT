#include<stdint.h>
#include <stdio.h>


#include<immintrin.h>

typedef uint64_t integer;
typedef __m256i integer_packed;



__m256i mul_mod_mersenne_avx256_64(__m256i a, __m256i b, __m256i p)
{
  /*

    I need to work with unsigned values, as there is no right
    arithmetic shift (with sign extension) for quad word - 64 bits integers.
    
    Each value occupies 64 bits (4 values per register) so that we
    can apply lazy reduction during split radix level2/3/4, by not
    reducing after add/sub and maybe also after multilication by
    sqrt(i) (left shift by 16).
    
    For multiplication, we need to reduce to 32 bits so we can use
    _mm256_mul_epu32 that multiplies 32bit integers at positions
    0,2,4,6 and save 64 bit results.
    
    
    Input: a and b should be <= 2^32-1. Do a double reduction before calling
    this subroutine
    
    
    Output: after 64 bit multiplication and one reduction step, the
    result is <= 10737418234 <= 5*p. This should be considered for
    subsequent lazy reductions.
    
  */    


  
  __m256i T, t1, t2, t3;

  T =  _mm256_mul_epu32(a, b);


  t1 = _mm256_and_si256(T, p);
 
  t2 = _mm256_srli_epi64(T, 31);
  
  t3 = _mm256_add_epi64(t1, t2);
  
  
  return t3;

}


__m256i mul_mod_mersenne_avx256_64_2(__m256i a, __m256i b, __m256i p)
{

  /*

     Input: a and b should be <= 2^32-1. Do a reduction before calling
     this subroutine


     Output: after 64 bit multiplication and two reduction steps, the result is in [0:p+5]

  */    

  

  __m256i T, t1, t2, t3;

   T =  _mm256_mul_epu32(a, b);

  t1 = _mm256_and_si256(T, p);
 
  t2 = _mm256_srli_epi64(T, 31);
  
  T = _mm256_add_epi64(t1, t2);

  t1 = _mm256_and_si256(T, p);
 
  t2 = _mm256_srli_epi64(T, 31);
  
  t3 = _mm256_add_epi64(t1, t2);
  
  return t3;

}


__m256i red_mersenne_avx256_64(__m256i T, __m256i p)
{

  __m256i t1, t2, t3;

  
  t1 = _mm256_and_si256(T, p);
 
  t2 = _mm256_srli_epi64(T, 31);
  
  t3 = _mm256_add_epi64(t1, t2);
    
  return t3;
  
}



__m256i red_mersenne_avx256_64_2(__m256i T, __m256i p)
{

  __m256i t1, t2, t3;

    
  
  t1 = _mm256_and_si256(T, p);
 
  t2 = _mm256_srli_epi64(T, 31);
  
  T = _mm256_add_epi64(t1, t2);

    
  t1 = _mm256_and_si256(T, p);
 
  t2 = _mm256_srli_epi64(T, 31);
  
  t3 = _mm256_add_epi64(t1, t2);

  return t3;
  
}



__m256i red_mersenne_avx256_64_2_final(__m256i T, __m256i p)
{

  __m256i t1, t2, t3, k;

    
  t1 = _mm256_and_si256(T, p);

  t2 = _mm256_srli_epi64(T, 31);

  T = _mm256_add_epi64(t1, t2);

  t1 = _mm256_and_si256(T, p);

  t2 = _mm256_srli_epi64(T, 31);

  t3 = _mm256_add_epi64(t1, t2);//x

  t1 = _mm256_sub_epi64(t3, p); //x-p


  //k = true if p > x)
  k = _mm256_cmpgt_epi64 (p, t3); //p > x

  //if k true select x, else x-p
  t2 =  _mm256_blendv_epi8(t1, t3, k);
  
  
  return t2;
  
}
