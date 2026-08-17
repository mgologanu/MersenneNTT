#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"


void finalize_r8_from_c4(uint64_t * x0,  uint64_t * y0, uint64_t *  x1, uint64_t * y1, uint64_t *  x2, uint64_t *   y2, uint64_t *   x3, uint64_t *  y3)
{
  uint64_t t0,t1,t2,t3,t4,t5;
  
  t0 = *x0 + *y0;
  t1 = *x0 - *y0;

  red1_scalar(&t0);
  *x0 = t0;

  red1_scalar(&t1);
  *y0 = t1;

  t0 = PRIME - *y2;
  red1_scalar(&t0);
  *y2 = t0;


  //  t0 = 0.5*(*x1 + *x3);
  t0 = *x1 + *x3 ;
  red1_scalar(&t0);
  t0 = t0 << 30;
  red1_scalar(&t0);
    
  // t1 = HALF_SQRT_HALF*(*x3 - *x1);
  t1 = *x3 - *x1;
  red1_scalar(&t1);
  t1 = t1 << 14;
  red1_scalar(&t1);


  //t2 = HALF_SQRT_HALF*(*y1 + *y3);
  t2 = *y1 + *y3;
  red1_scalar(&t2);
  t2 = t2 << 14;
  red1_scalar(&t2);
    
  //t3 = 0.5*(*y1 - *y3);
  t3 = *y1 - *y3;
  red1_scalar(&t3);
  t3 = t3 << 30;
  red1_scalar(&t3);
  

  //  A = t0 + i*t3
  //  B = t2 + i*t1

  //  D = B * (1-i)

  t4 = t1 + t2;
  t5 = t1 - t2;

  //A = t0 + i*t3
  //D = t4 + i*t5

  //z1 = A + D;

  //  *x1 = t0 + t4;
  t1 = t0 + t4;
  red1_scalar(&t1);
  *x1 = t1;

  //  *y1 = t3 + t5;
  t1 = t3 + t5;
  red1_scalar(&t1);
  *y1 = t1;

  //z3 = conj(A) - conj(D)
  //  *x3 = t0 - t4;
  t1 = t0 - t4;
  red1_scalar(&t1);
  *x3 = t1;
  
  //  *y3 = t5 - t3;
  t1 = t5 - t3;
  red1_scalar(&t1);
  *y3 = t1;
  
}

/* #define P1_F 0x20 */
/* #define P2_F 0x31 */

/* #define I1_F 0x0 */
/* #define I2_F 0xF */



void r32(register integer *a)
{
  integer_packed t0, t1, t2, t3;
  
  rpass_srl3(a, d16, d32, 4);//R(x^8-1) and 3 x C(x^4-1)

  /* 
     transform R8 = R(x^8-1) in C(x^4-1) for even and odd components

     R8 -> C4(even, odd)

     Need to shuffle:

     even(R8) -> integerC4) odd(R8) -> imag(C4)
  */
  
  t0 = LOAD(a[0]);
  t1 = LOAD(a[4]);

  // t2 = PERM(t0,t1,P1_F);
  // t3 = PERM(t0,t1,P2_F);

  t2 = _mm256_permute2f128_si256(t0,t1,0x20) ;
  t3 = _mm256_permute2f128_si256(t0,t1,0x31) ;

  // t0 = SHUF(t2,t3,I1_F);
  // t1 = SHUF(t2,t3,I2_F);

  
  t0 =  (__m256i ) _mm256_shuffle_pd((__m256d) t2,(__m256d) t3, 0x0);
  t1 =  (__m256i ) _mm256_shuffle_pd((__m256d) t2,(__m256d) t3, 0xF);
  
  STORE0(a[0], t0);
  STORE0(a[4], t1);

  /*
    Parallel (vertical) 4 x C4
    Note that c4_vert also contains the transpose part
  */
  
  c4_vert(a);

  /*
    Still not done - R8 from C4(even, odd) needs supplementary massaging
    
    Note that the output of c4_vert is ordered as {0, 2, 3, 1} (see freq.c)

    Also the 4xC4 is transposed, so that the frequencies {f0, f1, f2,
    f3} for the first C4 are in a[0], a[24], a[8], a[16].

  */
  
   finalize_r8_from_c4(&a[0],&a[4],&a[24],&a[28],&a[8],&a[12],&a[16],&a[20]);
  
}

void  r64(register integer *a)
{
  /* 
     Could do SRL3 and obtain R16 and 3xC8, then R16->C8, 4xC8 and finalize for R16 from C8(even, odd).
     
     A simpler solution:

     SR -> R32 and 1xC16

     Even if this is slower than SRL3, it is applied only once in r2^n, so it
     adds a constant time for all values of n.
  */
  rpass_sr(a, d64, 8);
  r32(a);
  cs16(a+32);
}

void r128(register integer *a)
{
  rpass_srl3(a, d64, d128, 16);

  r32(a);
  cs16(a+32);
  cs16(a+64);
  cs16(a+96);
}


void r256(register integer *a)
{
  rpass_srl3(a, d128, d256, 32);
  
  r64(a);
  cs32(a + 64);
  cs32(a + 128);
  cs32(a + 192);
}


