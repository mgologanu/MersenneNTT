
#include "ntt.h"
#include "pre8_V4.h"

#include "i.h"

void transpose(integer_packed *x0, integer_packed *x1, integer_packed *x2, integer_packed *x3, integer * a) 

{
  integer_packed t0, t1, t2, t3;

  size_t offset = 0;
  size_t stride = 8;

  t0 = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_load_si128((__m128i *) &a[offset+0*stride+0])), _mm_load_si128((__m128i *) &a[offset+2*stride+0]), 1); 
  t1 = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_load_si128((__m128i *) &a[offset+1*stride+0])), _mm_load_si128((__m128i *) &a[offset+3*stride+0]), 1); 
  t2 = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_load_si128((__m128i *) &a[offset+0*stride+2])), _mm_load_si128((__m128i *) &a[offset+2*stride+2]), 1); 
  t3 = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_load_si128((__m128i *) &a[offset+1*stride+2])), _mm_load_si128((__m128i *) &a[offset+3*stride+2]), 1); 
  
  *x0 = _mm256_unpacklo_epi64(t0,t1);					
  *x1 = _mm256_unpackhi_epi64(t0,t1);		
  *x2 = _mm256_unpacklo_epi64(t2,t3);		
  *x3 = _mm256_unpackhi_epi64(t2,t3);

  return;
}
