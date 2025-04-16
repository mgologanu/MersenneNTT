#include "ntt.h"
#include "pre8_V4.h"

#include "i.h"

/* 
   FFT4V - 4*V complex numbers

   V <= 4

   a: 4*V complex = 8*V values  interleaved V-packed format  - V real parts, V imag parts, etc.


   transpose and apply FFT4 using vertical vectorization

   Input:
    -C4-     -C4-    -C4-   -C4-
   x0 y0   x1 y1   x2 y2   x3 y3 
    0  V   2V 3V   4V 5V   6V 7V 

    - transpose(x0, x1, x2, x3)

    - transpose(y0, y1, y2, y3)

    - Apply C4 to x and y


    Output order: 0 2 3 1, see freq.c

    Note: the partially split format, given by the output of rad4 or
    rad8 (twice) is precisely the same - there is no need for a cs4_vert function.

*/
void c4_vert(integer *a)
{
  /* const integer_packed  p = BROADCAST(PRIME); */

  /* const integer_packed p_64 = BROADCAST(PRIME_64); */

  /* integer_packed _red1_t1, _red1_t2, _red1_t3; */
  /* integer_packed _t0, _t1, _t2, _t3, retval; */
  /* integer_packed  _AipBr, _AimBr, _ArpBi, _ArmBi; */
  /* integer_packed _Ar, _Br, _Ai, _Bi, _t02, _t13, _s02, _s13; */
  
  integer_packed x0, y0, x1, y1, x2, y2, x3, y3;

  /* integer_packed _s0, _s1, _s2, _s3; */

  /* x0 = LOAD(a[0*V]); */
  /* y0 = LOAD(a[1*V]); */

  /* x1 = LOAD(a[2*V]); */
  /* y1 = LOAD(a[3*V]); */

  /* x2 = LOAD(a[4*V]); */
  /* y2 = LOAD(a[5*V]); */

  /* x3 = LOAD(a[6*V]); */
  /* y3 = LOAD(a[7*V]); */


  /* TRANSPOSE4(x0,x1,x2,x3); */
  /* TRANSPOSE4(y0,y1,y2,y3); */


  /* TRANSPOSE4_FROM_MEM(x0,x1,x2,x3,a,0,V2); */
  /* TRANSPOSE4_FROM_MEM(y0,y1,y2,y3,a,V,V2); */


  /* C4_REG_MEM(x0, y0, x1, y1, x2, y2, x3, y3, a[0*V], a[1*V], a[2*V], a[3*V], a[4*V], a[5*V], a[6*V], a[7*V]); */
  transpose4_from_mem(&x0,&x1,&x2,&x3,a,0,V2);
  transpose4_from_mem(&y0,&y1,&y2,&y3,a,V,V2);

  
  c4_reg_to_mem(x0, y0, x1, y1, x2, y2, x3, y3, a, a+1*V, a+2*V, a+3*V, a+4*V, a+5*V, a+6*V, a+7*V);
  
}

