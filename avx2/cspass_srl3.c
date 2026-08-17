#include "ntt.h"

#include "pre8_V4.h"
#include "i.h"





/* 
   Extended split radix (or level 3) pass for x^8n-1 -> 8n complex numbers 

   a: 8n complex = 16n values = split real/imag arrays. Each array
   is contiguous and imag_array = real_array + 8n (full array is also
   contiguous)

   w: n complex numbers or 2n values = [0...2n-1] in interleaved V-packed format - 1/4th of roots of order 4n

   w2: n complex numbers or 2n values = [0...2n-1] in interleaved V-packed format - 1/8th of roots of order 8n

   No mirroring for w or w2!

*/

void cspass_srl3(integer *a,  const integer *w,  const integer *w2, size_t n)
{

  // const integer_packed  p = BROADCAST(PRIME);

  const integer_packed p_64 = BROADCAST(PRIME_64);
  

  /* integer_packed _mul1_t1, _mul1_t2, _mul1_t3, retval; */
  /* integer_packed _mul0_t1, _mul0_t2, _mul0_t3; */
  /* integer_packed _red1_t1, _red1_t2, _red1_t3; */
  
  /* integer_packed _t0, _t1, _t2,_AipBr, _AimBr, _ArpBi, _ArmBi; */
    
  integer_packed x0, x1, x2, x3, x4, x5, x6, x7;
  integer_packed y0, y1, y2, y3, y4, y5, y6, y7;
  integer_packed xp04, xp26, xp15, xp37, yp04, yp26, yp15, yp37; 
  integer_packed c, s, c2, s2, c3, s3;
  integer_packed Ar, Br, Ai, Bi;
  integer_packed t0, t1, t2, t3;
  integer_packed Er, Ei, Fr, Fi;
  integer_packed ar,br,cr,dr,ai,bi,ci,di;

  size_t k;

   
  integer *a1 = a  + n;
  integer *a2 = a1 + n;
  integer *a3 = a2 + n;
  integer *a4 = a3 + n;
  integer *a5 = a4 + n;
  integer *a6 = a5 + n;
  integer *a7 = a6 + n;
  
  integer *b  = a  + (n<<3);
  
  integer *b1 = b  + n;
  integer *b2 = b1 + n;
  integer *b3 = b2 + n;
  integer *b4 = b3 + n;
  integer *b5 = b4 + n;
  integer *b6 = b5 + n;
  integer *b7 = b6 + n;
  
  k = n >> VPOW;  //n/V

  do {
    
    x0 = LOAD(a[0]);
    x1 = LOAD(a1[0]);
    x2 = LOAD(a2[0]);
    x3 = LOAD(a3[0]);
    x4 = LOAD(a4[0]);
    x5 = LOAD(a5[0]);
    x6 = LOAD(a6[0]);
    x7 = LOAD(a7[0]);
    
    y0 = LOAD(b[0]);
    y1 = LOAD(b1[0]);
    y2 = LOAD(b2[0]);
    y3 = LOAD(b3[0]);
    y4 = LOAD(b4[0]);
    y5 = LOAD(b5[0]);
    y6 = LOAD(b6[0]);
    y7 = LOAD(b7[0]);
    
    xp04 = ADD(x0, x4);
    ar   = SUB(x0, x4);
    
    xp26 = ADD(x2, x6);
    br   = SUB(x2, x6);

    xp15 = ADD(x1, x5);
    cr   = SUB(x1, x5);
    
    xp37 = ADD(x3, x7);
    dr   = SUB(x3, x7);

    t0 = ADD(xp04, xp26); //X 0+4+2+6 = X0
    Ar = SUB(xp04, xp26); //X 0+4-2-6

    STORE(a[0], t0);

    t0 = ADD(xp15, xp37); //X 1+5+3+7 = X1
    Br = SUB(xp15, xp37); //X 1+5-3-7

    STORE(a1[0], t0);
    
    yp04 = ADD(y0, y4);
    ai   = SUB(y0, y4);
    
    yp26 = ADD(y2, y6);
    bi   = SUB(y2, y6);

    yp15 = ADD(y1, y5);
    ci   = SUB(y1, y5);
    
    yp37 = ADD(y3, y7);
    di   = SUB(y3, y7);

    t0 = ADD(yp04, yp26); //Y 0+4+2+6 = Y0
    Ai = SUB(yp04, yp26); //Y 0+4-2-6
    
    STORE(a2[0], t0);

    t0 = ADD(yp15, yp37); //Y 1+5+3+7 = Y1
    Bi = SUB(yp15, yp37); //Y 1+5-3-7
    
    STORE(a3[0], t0);

    c = LOAD(w[0]);
    s = LOAD(w[V]);

    //    TWISTED(Ar,Ai,Br,Bi,c,s,a4[0],a5[0],a6[0],a7[0]);

    twisted(Ar, Ai, Br, Bi, c, s, a4, a5, a6, a7);
 

    //Multiplication by 2^15, as sqrt(i) = 2^15*(1+i)
    cr = ADD(cr, p_64);
    cr = _mm256_slli_epi64(cr, SQ2_POW);
    red1(&cr);

    dr = ADD(dr, p_64);
    dr = _mm256_slli_epi64(dr, SQ2_POW);
    red1(&dr);

    ci = ADD(ci, p_64);
    ci = _mm256_slli_epi64(ci, SQ2_POW);
    red1(&ci);

    di = ADD(di, p_64);
    di = _mm256_slli_epi64(di, SQ2_POW);
    red1(&di);


    t0 = ADD(cr, dr); //cr + dr
    t1 = SUB(cr, dr); //cr - dr

    Ar = ADD(ar, t1); //ar +  cr - dr
    Er = SUB(ar, t1); //ar - (cr - dr)
    
    Br = ADD(br, t0); // br + cr + dr
    Fr = SUB(t0, br); //-br + cr + dr
    
    t0 = ADD(ci, di); //ci + di
    t1 = SUB(ci, di); //ci - di

    Ai = ADD(ai, t1); //ai +  ci - di
    Ei = SUB(ai, t1); //ai - (ci - di)

    Bi = ADD(bi, t0); // bi + ci + di
    Fi = SUB(t0, bi); //-bi + ci + di

    c2 = LOAD(w2[0]);
    s2 = LOAD(w2[V]);

    //    TWISTED(Ar,Ai,Br,Bi,c2,s2, b[0],b1[0],b6[0],b7[0]);

    twisted(Ar, Ai, Br, Bi, c, s, b, b1, b6, b7);

    // c3 + i s3 = (c + i s) (c2 + i s2)
    // The roots are already reduced and positive
    t0 = mul0(c, c2);
    t1 = mul0(s, s2);
    t2 = mul0(c, s2);
    t3 = mul0(s, c2);

    c3 = SUB(t0, t1);
    s3 = ADD(t2, t3);

    red1(&c3);
    red1(&s3);

    //   TWISTED(Er,Ei,Fr,Fi,c3,s3, b4[0],b5[0],b2[0],b3[0]);

    twisted(Er, Ei, Fr, Fi, c, s, b4, b5, b2, b3);

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
    
  }  while (k -= 1);
}

