#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

//#include <stdio.h>
void scalec(integer *a, size_t n)
{
  integer_packed pa;
  int n_pow = -1;
  size_t nn = n;

  
  while (nn>0) {
    nn = nn>>1;
    n_pow++;
 
  };

  n_pow = P_POW - n_pow;


  for (int i = 0; i<(n>>(VPOW-1)); i++)
    {
      pa = _mm256_load_si256((__m256i *) a);
      red1(&pa);
      pa = _mm256_slli_epi64(pa, n_pow);
      red1(&pa);
      _mm256_store_si256((__m256i *) a, pa);
      a += V;
    }
  

  
}
