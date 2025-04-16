#include "ntt.h"
#include "pre8_V4.h"

#include "i.h"



/* 

   Inverse transform for ts16_V4

   Recall t16_V4:

   Transpose 16 x V complex numbers to V x 16 appropriate as input for
   vertical C16

   Partially split arrays for real and imaginary parts as given by cspass_rad4 or cspass_rad8 (which actually gives 8xC16)

   --------- C16 ----------  --------- C16 ----------  ------------- C16 ----------  ------------- C16 -------------- 
   x0 x1 x2 x3  y0 y1 y2 y3  x4 x5 x6 x7  y4 y5 y6 y7  x8 x9 x10 x11  y8 y9 y10 y11  x12 x13 x14 x15  y12 y13 y14 y15 
   
    0  1  2  3  4  5  6  7   8  9 10  11  12 13 14 15  16 17 18  19  20 21  22  23   24  25  26  27   28  29  30  31   

   Each x and y is a packed array with V values

   transpose(x0 x4  x8 x12) -> save in x0   x1  x2  x3
   transpose(x1 x5  x9 x13) -> save in x4   x5  x6  x7
   transpose(x2 x6 x10 x14) -> save in x8   x9 x10 x11
   transpose(x3 x7 x11 x15) -> save in x12 x13 x14 x15

   Repeat for y

   INVERSE:

   transpose( x0  x1  x2  x3)  -> save in x0 x4  x8 x12
   transpose( x4  x5  x6  x7)  -> save in x1 x5  x9 x13
   transpose( x8  x9 x10 x11)  -> save in x2 x6 x10 x14
   transpose(x12 x13 x14 x15)  -> save in x3 x7 x11 x15
   
   Repeat for y

*/



void uts16_V4(integer *a)
{
  integer_packed t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15;
  // integer_packed s0, s1, s2, s3;

  /* t0 =  LOAD(a[0]); */
  /* t1 =  LOAD(a[V]); */
  /* t2 =  LOAD(a[2*V]); */
  /* t3 =  LOAD(a[3*V]); */

  /* TRANSPOSE4(t0,t1,t2,t3); */

  transpose4_from_mem(&t0,&t1,&t2,&t3, a, 0, V);
   
  STORE(a[0],    t0);

  /* t4 =  LOAD(a[8*V]); */
  /* t5 =  LOAD(a[9*V]); */
  /* t6 =  LOAD(a[10*V]); */
  /* t7 =  LOAD(a[11*V]); */

  /* TRANSPOSE4(t4,t5,t6,t7); */

  transpose4_from_mem(&t4,&t5,&t6,&t7, a, 8*V, V);

  STORE(a[8*V],  t1);
  
  /* t8 =  LOAD(a[16*V]); */
  /* t9 =  LOAD(a[17*V]); */
  /* t10 =  LOAD(a[18*V]); */
  /* t11 =  LOAD(a[19*V]); */

  /* TRANSPOSE4(t8,t9,t10,t11); */

  transpose4_from_mem(&t8,&t9,&t10,&t11, a, 16*V, V);

  STORE(a[16*V],  t2);

  /* t12 =  LOAD(a[24*V]); */
  /* t13 =  LOAD(a[25*V]); */
  /* t14 =  LOAD(a[26*V]); */
  /* t15 =  LOAD(a[27*V]); */

  /* TRANSPOSE4(t12,t13,t14,t15); */


  transpose4_from_mem(&t12,&t13,&t14,&t15, a, 24*V, V);
   
  STORE0(a[24*V], t3);
  STORE0(a[1*V],  t4);
  STORE0(a[9*V],  t5);
  STORE0(a[17*V],  t6);
  STORE0(a[25*V], t7);
  STORE0(a[2*V],  t8);
  STORE0(a[10*V],  t9);
  STORE0(a[18*V], t10);
  STORE0(a[26*V], t11);
  STORE0(a[3*V],  t12);
  STORE0(a[11*V],  t13);
  STORE0(a[19*V], t14);
  STORE0(a[27*V], t15);

  //Same for y
  a += V << 2;

  /* t0 =  LOAD(a[0]); */
  /* t1 =  LOAD(a[V]); */
  /* t2 =  LOAD(a[2*V]); */
  /* t3 =  LOAD(a[3*V]); */

  /* TRANSPOSE4(t0,t1,t2,t3); */

  transpose4_from_mem(&t0,&t1,&t2,&t3, a, 0, V);

  STORE(a[0],    t0);

  /* t4 =  LOAD(a[8*V]); */
  /* t5 =  LOAD(a[9*V]); */
  /* t6 =  LOAD(a[10*V]); */
  /* t7 =  LOAD(a[11*V]); */

  /* TRANSPOSE4(t4,t5,t6,t7); */

  transpose4_from_mem(&t4,&t5,&t6,&t7, a, 8*V, V);
    
  STORE(a[8*V],  t1);
  
  /* t8 =  LOAD(a[16*V]); */
  /* t9 =  LOAD(a[17*V]); */
  /* t10 =  LOAD(a[18*V]); */
  /* t11 =  LOAD(a[19*V]); */

  /* TRANSPOSE4(t8,t9,t10,t11); */

  transpose4_from_mem(&t8,&t9,&t10,&t11, a, 16*V, V);
  
  STORE(a[16*V],  t2);

  /* t12 =  LOAD(a[24*V]); */
  /* t13 =  LOAD(a[25*V]); */
  /* t14 =  LOAD(a[26*V]); */
  /* t15 =  LOAD(a[27*V]); */

  /* TRANSPOSE4(t12,t13,t14,t15); */

  transpose4_from_mem(&t12,&t13,&t14,&t15, a, 24*V, V);
 
  STORE0(a[24*V], t3);
  STORE0(a[1*V],  t4);
  STORE0(a[9*V],  t5);
  STORE0(a[17*V],  t6);
  STORE0(a[25*V], t7);
  STORE0(a[2*V],  t8);
  STORE0(a[10*V],  t9);
  STORE0(a[18*V], t10);
  STORE0(a[26*V], t11);
  STORE0(a[3*V],  t12);
  STORE0(a[11*V],  t13);
  STORE0(a[19*V], t14);
  STORE0(a[27*V], t15);
  
}
