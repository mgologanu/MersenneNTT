#include<stdint.h>
#include <stdio.h>


#include<immintrin.h>

typedef int64_t integer;
typedef __m512i integer_packed;



integer_packed mul_mod_mersenne_avx512_64(integer_packed a, integer_packed b, integer_packed p)
{

  // assuming that each 64 bit integer actually occupies only 32 bit max.
  //Note that the sign bit is extended up to 64 bit

  integer_packed T, t1, t2, t3;

  T = _mm512_mullo_epi64(a, b);
  
  t1 = _mm512_and_epi64(T, p);

  t2 = _mm512_srai_epi64(T, 31);

  t3 = _mm512_add_epi64(t1, t2);
  
  return t3;

}


integer_packed mul_mod_mersenne_avx512_64_2(integer_packed a, integer_packed b, integer_packed p)
{

  // assuming that each 64 bit integer actually occupies only 32 bit max.
  //Note that the sign bit is extended up to 64 bit

  integer_packed T, t1, t2, t3;

  T = _mm512_mullo_epi64(a, b);
  
  t1 = _mm512_and_epi64(T, p);

  t2 = _mm512_srai_epi64(T, 31);

  T = _mm512_add_epi64(t1, t2);

  t1 = _mm512_and_epi64(T, p);

  t2 = _mm512_srai_epi64(T, 31);

  t3 = _mm512_add_epi64(t1, t2);
    
  return t3;

}


integer_packed red_mersenne_avx512_64(integer_packed T, integer_packed p)
{

  integer_packed t1, t2, t3;

  t1 = _mm512_and_epi64(T, p);

  t2 = _mm512_srai_epi64(T, 31);

  t3 = _mm512_add_epi64(t1, t2);
    
  return t3;
  
}



integer_packed red_mersenne_avx512_64_2(integer_packed T, integer_packed p)
{

  integer_packed t1, t2, t3;

  t1 = _mm512_and_epi64(T, p);

  t2 = _mm512_srai_epi64(T, 31);

  T = _mm512_add_epi64(t1, t2);

  t1 = _mm512_and_epi64(T, p);

  t2 = _mm512_srai_epi64(T, 31);

  t3 = _mm512_add_epi64(t1, t2);
    
  return t3;
  
}



integer_packed red_mersenne_avx512_64_2_final(integer_packed T, integer_packed p, integer_packed pm1_half)
{

  integer_packed t1, t2, t3;

  __mmask8 k;
    
  t1 = _mm512_and_epi64(T, p);

  t2 = _mm512_srai_epi64(T, 31);

  T = _mm512_add_epi64(t1, t2);

  t1 = _mm512_and_epi64(T, p);

  t2 = _mm512_srai_epi64(T, 31);

  t3 = _mm512_add_epi64(t1, t2);

  t1 = _mm512_sub_epi64(t3, p);


  //k = true if pm1_half < t3)
  k = _mm512_cmp_epi64_mask (pm1_half, t3, _MM_CMPINT_LT);

  //if k true select t1, else t3
  t2 =  _mm512_mask_blend_epi64 (k, t3, t1);
  
  
  return t2;
  
}
