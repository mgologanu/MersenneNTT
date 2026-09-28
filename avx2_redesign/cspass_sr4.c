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


void cspass_sr4(int32_t *a,  const int32_t  *w, size_t n)
{

  __m256i x0, x1, x2, x3, y0, y1, y2, y3;
  __m256i c, s;
  //  __m256i Ar, Br, Ai, Bi;
  // __m256i t02, t13, s02, s13;

  size_t k, n2;

  n2 = n << 1; //2n
  
  int32_t *a1 = a  + n2;
  int32_t *a2 = a1 + n2;
  int32_t *a3 = a2 + n2;
  
  int32_t *b  = a  + (n<<3);
  
  int32_t *b1 = b  + n2;
  int32_t *b2 = b1 + n2;
  int32_t *b3 = b2 + n2;

  
  k = n2 >> 3;  //n/V where V = avx2_vector_length

  
  do    {

    x0 = _mm256_load_si256((__m256i *) &a[0]);
    x1 = _mm256_load_si256((__m256i *) &a1[0]);
    x2 = _mm256_load_si256((__m256i *) &a2[0]);
    x3 = _mm256_load_si256((__m256i *) &a3[0]);

    
    y0 = _mm256_load_si256((__m256i *) &b[0]);
    y1 = _mm256_load_si256((__m256i *) &b1[0]);
    y2 = _mm256_load_si256((__m256i *) &b2[0]);
    y3 = _mm256_load_si256((__m256i *) &b3[0]);


    c = _mm256_load_si256((__m256i *) &w[0]);
    s = _mm256_load_si256((__m256i *) &w[V]);


    add_sub(&x0, &x2);
    add_sub(&x1, &x3);

    add_sub(&y0, &y2);
    add_sub(&y1, &y3);
       

    // t02 = add(x0, x2); //X0
    //    Ar  = sub(x0, x2);
    
    // t13 = add(x1, x3); //X1
    // Br  = sub(x1, x3);
    
    // s02 = add(y0, y2); //Y0
    // Ai  = sub(y0, y2);
    
    // s13 = add(y1, y3); //Y1
    // Bi  = sub(y1, y3);

    _mm256_store_si256((__m256i *) &a[0],  x0);

    _mm256_store_si256((__m256i *) &a1[0], x1);

    _mm256_store_si256((__m256i *) &a2[0], y0);

    _mm256_store_si256((__m256i *) &a3[0], y1);
    
    //apply twister and save to array
    
    //   twister(Ar, Ai, Br, Bi, c, s, b, b1, b2, b3);

    twister(x2, y2, x3, y3, c, s, b, b1, b2, b3);
      
    a  += V;
    a1 += V;
    a2 += V;
    a3 += V;
 
    b  += V;
    b1 += V;
    b2 += V;
    b3 += V;

    w  += V2;

    
  }  while (k -= 1);
}
