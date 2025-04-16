#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"


void mulr(integer *a, integer *b, size_t n)
{

  integer t0, t1;
  //save the real values.
  t0 = a[0] * b[0];
  t1 = a[V] * b[V];
  
  mulc(a, b, n>>1);

  //put the real values back
  a[0] = t0;
  a[V] = t1;
  
}
