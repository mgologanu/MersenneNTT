#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"


void mulc(register integer *a, register integer *b, size_t n)
{
  integer_packed t1, t2, t3, t4, t5, t6;

  size_t k;
  
  for (k = 0; k < (n>>VPOW); k++)
    {
      t1 = LOAD(a[0]);
      t2 = LOAD(a[V]);
      
      t3 = LOAD(b[0]);
      t4 = LOAD(b[V]);
      
      //mul1 requires second argument less than 2*p
      red1(&t3);
      red1(&t4);
      
      t5 = mul1(t1, t3) - mul1(t2, t4);
      t6 =  mul1(t1, t4) + mul1(t2, t3);
      
      STORE(a[0], t5);
      STORE(a[V], t6);
      
      a += V2;
      b += V2;
  }
  
}
