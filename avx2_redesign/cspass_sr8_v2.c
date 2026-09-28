#include <stdlib.h>
#include <stdint.h>
#include<immintrin.h>


#include "mrsn_avx2_internal.h"

#define V 8
#define V2 16

/* 
   Extended split radix (or level 3) pass for 8n complex values, or an
   array of lengh 2*8*n = 16*n, with split real and imaginary
   values (first half real values, second half imaginary values).

   a: 8n complex = 16n values = split real/imag arrays.

   w: n complex numbers or  2n array = [0...2n-1] in interleaved V-packed format  - 1/4th of roots of order 4n

   w2: n complex numbers or 2n values = [0...2n-1] in interleaved V-packed format - 1/8th of roots of order 8n

   No mirroring for w or w2!

*/


void cspass_sr8(int32_t *a,  const int32_t  *w,  const int32_t *w2, const int32_t *w3, size_t n)
{

  __m256i x0, x1, x2, x3, x4, x5, x6, x7;
  __m256i y0, y1, y2, y3, y4, y5, y6, y7;
  __m256i xp04, xp26, xp15, xp37, yp04, yp26, yp15, yp37; 
  __m256i c, s, c2, s2, c3, s3;
  __m256i Ar, Br, Ai, Bi;
  __m256i t0, t1, t2, t3;
  __m256i Er, Ei, Fr, Fi;
  __m256i ar,br,cr,dr,ai,bi,ci,di;

  size_t k;
  
  int32_t *a1 = a  + n;
  int32_t *a2 = a1 + n;
  int32_t *a3 = a2 + n;
  int32_t *a4 = a3 + n;
  int32_t *a5 = a4 + n;
  int32_t *a6 = a5 + n;
  int32_t *a7 = a6 + n;
  
  int32_t *b  = a  + (n<<3);
  
  int32_t *b1 = b  + n;
  int32_t *b2 = b1 + n;
  int32_t *b3 = b2 + n;
  int32_t *b4 = b3 + n;
  int32_t *b5 = b4 + n;
  int32_t *b6 = b5 + n;
  int32_t *b7 = b6 + n;
  
  k = n >> 3;  //n/V where V = avx2_vector_length

  
  do    {

    x0 = _mm256_load_si256((__m256i *) &a[0]);
    x1 = _mm256_load_si256((__m256i *) &a1[0]);
    x2 = _mm256_load_si256((__m256i *) &a2[0]);
    x3 = _mm256_load_si256((__m256i *) &a3[0]);
    x4 = _mm256_load_si256((__m256i *) &a4[0]);
    x5 = _mm256_load_si256((__m256i *) &a5[0]);
    x6 = _mm256_load_si256((__m256i *) &a6[0]);
    x7 = _mm256_load_si256((__m256i *) &a7[0]);
    
    y0 = _mm256_load_si256((__m256i *) &b[0]);
    y1 = _mm256_load_si256((__m256i *) &b1[0]);
    y2 = _mm256_load_si256((__m256i *) &b2[0]);
    y3 = _mm256_load_si256((__m256i *) &b3[0]);
    y4 = _mm256_load_si256((__m256i *) &b4[0]);
    y5 = _mm256_load_si256((__m256i *) &b5[0]);
    y6 = _mm256_load_si256((__m256i *) &b6[0]);
    y7 = _mm256_load_si256((__m256i *) &b7[0]);

    
    xp04 = add(x0, x4);
    ar   = sub(x0, x4);
  
    xp26 = add(x2, x6);
    br   = sub(x2, x6);
  
    xp15 = add(x1, x5);
    cr   = sub(x1, x5);
  
    xp37 = add(x3, x7);
    dr   = sub(x3, x7);

    t0 = add(xp04, xp26); //X 0+4+2+6 = X0
    Ar = sub(xp04, xp26); //X 0+4-2-6
  
    _mm256_store_si256((__m256i *) &a[0], t0);
  
    t0 = add(xp15, xp37); //X 1+5+3+7 = X1
    Br = sub(xp15, xp37); //X 1+5-3-7
  
    _mm256_store_si256((__m256i *) &a1[0], t0);

    yp04 = add(y0, y4);
    ai   = sub(y0, y4);
    
    yp26 = add(y2, y6);
    bi   = sub(y2, y6);
  
    yp15 = add(y1, y5);
    ci   = sub(y1, y5);
  
    yp37 = add(y3, y7);
    di   = sub(y3, y7);
  
    t0 = add(yp04, yp26); //Y 0+4+2+6 = Y0
    Ai = sub(yp04, yp26); //Y 0+4-2-6
  
    _mm256_store_si256((__m256i *) &a2[0], t0);

    t0 = add(yp15, yp37); //Y 1+5+3+7 = Y1
    Bi = sub(yp15, yp37); //Y 1+5-3-7

    _mm256_store_si256((__m256i *) &a3[0], t0);


    c = _mm256_load_si256((__m256i *) &w[0]);
    s = _mm256_load_si256((__m256i *) &w[V]);

    //apply twister and save to array
    
    twister2(Ar, Ai, Br, Bi, c, s, a4, a5, a6, a7);

    //Multiplication by 2^15, as sqrt(i) = 2^15*(1+i)

    cr = rot15(cr);
    dr = rot15(dr);
    ci = rot15(ci);
    di = rot15(di);

    t0 = add(cr, dr); //cr + dr
    t1 = sub(cr, dr); //cr - dr

    Ar = add(ar, t1); //ar +  cr - dr
    Er = sub(ar, t1); //ar - (cr - dr)
    
    Br = add(br, t0); // br + cr + dr
    Fr = sub(t0, br); //-br + cr + dr
    
    t0 = add(ci, di); //ci + di
    t1 = sub(ci, di); //ci - di

    Ai = add(ai, t1); //ai +  ci - di
    Ei = sub(ai, t1); //ai - (ci - di)

    Bi = add(bi, t0); // bi + ci + di
    Fi = sub(t0, bi); //-bi + ci + di

    c2 = _mm256_load_si256((__m256i *) &w2[0]);
    s2 = _mm256_load_si256((__m256i *) &w2[V]);

    twister2(Ar, Ai, Br, Bi, c2, s2, b, b1, b6, b7);

    /*
    t0 = mul(c, c2);
    t1 = mul(s, s2);
    t2 = mul(c, s2);
    t3 = mul(s, c2);
    
    c3 = sub(t0, t1);
    s3 = add(t2, t3);
    */

    
    c3 = _mm256_load_si256((__m256i *) &w3[0]);
    s3 = _mm256_load_si256((__m256i *) &w3[V]);

    twister2(Er, Ei, Fr, Fi, c3, s3, b4, b5, b2, b3);
      
    a  += V;
    a1 += V;
    a2 += V;
    a3 += V;
    a4 += V;
    a5 += V;
    a6 += V;
    a7 += V;

    b  += V;
    b1 += V;
    b2 += V;
    b3 += V;
    b4 += V;
    b5 += V;
    b6 += V;
    b7 += V;

    w  += V2;
    w2 += V2;
    w3 += V2;
    
  }  while (k -= 1);
}
