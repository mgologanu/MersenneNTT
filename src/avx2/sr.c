#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void scaler(integer *a, size_t n)
{

  integer_packed pa;

  integer s0,s1;

  integer *a1 = a;
  
  int n_pow = -1;
  size_t nn = n;
  
  while (nn>0) {
    nn = nn>>1;
    n_pow++;
  };

  n_pow = P_POW - n_pow;

  s0 = a[0];
  s1 = a[V];

  
  red1_scalar(&s0);
  s0 = s0 << n_pow;
  red1_scalar(&s0);

  red1_scalar(&s1);
  s1 = s1 << n_pow;
  red1_scalar(&s1);
  
  n_pow = n_pow+1;

  for (int i = 0; i<(n>>VPOW); i++)
    {
      pa = _mm256_load_si256((__m256i *) a);
      red1(&pa);
      pa = _mm256_slli_epi64(pa, n_pow);
      red1(&pa);
      _mm256_store_si256((__m256i *) a, pa);
      a += V;
    }


  a1[0] = s0;
  a1[V] = s1;

  

  
}
