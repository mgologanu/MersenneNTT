#include<stdint.h>
#include <stdio.h>


#include<immintrin.h>

typedef int64_t integer;
typedef __m256i integer_packed;



int64_t mul_mod_mersenne_avx512_64(integer_packed a, integer_packed b, integer_packed p)
{

  // assuming that each 64 bit integer actually occupies only 32 bit max.
  //Note that The sign bit is extended up to 64 bit

  integer_packed T, t1;

  T = _mm256_mul_epi32(a, b);
  
  t1 = _mm256_and_si256(T, p);

  

  t2 =  (T >> p_pow);

  t3 = t1 + t2;
    
  // printf("%ld %ld %ld %ld\n", T, t1, t2, t3);
  
  return t3;

}
