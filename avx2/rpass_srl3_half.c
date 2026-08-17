#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"




/* 
   Extended split radix (or level 3) pass for x^8n-1 - 8n real numbers 

   a: 8n real = 8n values in normal order (unpacked)

   w: n complex numbers or 2n values = [0...2n-1] in interleaved V-packed format - 1/4th of roots of order 4n

   w2: n complex numbers or 2n values = [0...2n-1] in interleaved V-packed format - 1/8th of roots of order 8n

   No mirroring for w!

*/
void rpass_srl3_half(integer *a, const integer *w, const integer *w2, size_t n)
{
  const integer_packed p_64 = BROADCAST(PRIME_64);

  integer_packed x0, x1, x2, x3; 
  integer_packed xp04, xp26, xp15, xp37;
  integer_packed c, s, c2, s2, c3, s3;
  integer_packed Ar, Br;
  integer_packed t0, t1, t2, t3;
  integer_packed Er, Fr;
  integer_packed ar,br,cr,dr;

  size_t k;
    
  integer *a1;
  integer *a2;
  integer *a3;
  integer *a4;
  integer *a5;
  integer *a6;
  integer *a7;
  
   
  a1 = a  + n;
  a2 = a1 + n;
  a3 = a2 + n;
  a4 = a3 + n;
  a5 = a4 + n;
  a6 = a5 + n;
  a7 = a6 + n;
  
  k = n >> VPOW;  //n/V

  do {
    
    x0 = LOAD(a[0]);
    x1 = LOAD(a1[0]);
    x2 = LOAD(a2[0]);
    x3 = LOAD(a3[0]);
    /* x4 = LOAD(a4[0]); */
    /* x5 = LOAD(a5[0]); */
    /* x6 = LOAD(a6[0]); */
    /* x7 = LOAD(a7[0]); */
    
    /* xp04 = ADD(x0, x4); */
    /* ar   = SUB(x0, x4); */
    
    /* xp26 = ADD(x2, x6); */
    /* br   = SUB(x2, x6); */

    /* xp15 = ADD(x1, x5); */
    /* cr   = SUB(x1, x5); */
    
    /* xp37 = ADD(x3, x7); */
    /* dr   = SUB(x3, x7); */

    
    xp04 = x0; 
    ar   = x0; 
    
    xp26 = x2; 
    br   = x2; 

    xp15 = x1; 
    cr   = x1; 
    
    xp37 = x3; 
    dr   = x3; 
  
      
    
    t0 = ADD(xp04, xp26); //X 0+4+2+6 = X0
    Ar = SUB(xp04, xp26); //X 0+4-2-6
    STORE(a[0], t0);

    t0 = ADD(xp15, xp37); //X 1+5+3+7 = X1
    Br = SUB(xp15, xp37); //X 1+5-3-7
    STORE(a1[0], t0);

    c = LOAD(w[0]);
    s = LOAD(w[V]);

    //    REAL_TWISTED(Ar,Br,c,s,a2[0],a3[0]);

    real_twisted(Ar, Br, c, s, a2, a3);
 
    /* q2 = BROADCAST(sq2[0]); */

    /* cr = MUL(q2, cr); */
    /* dr = MUL(q2, dr); */

    //Multiplication by 2^15, as sqrt(i) = 2^15*(1+i)
    cr = ADD(cr, p_64);
    cr = _mm256_slli_epi64(cr, SQ2_POW);
    red1(&cr);

    dr = ADD(dr, p_64);
    dr = _mm256_slli_epi64(dr, SQ2_POW);
    red1(&dr);


    
    
    t0 = ADD(cr, dr); //cr + dr
    t1 = SUB(cr, dr); //cr - dr

    Ar = ADD(ar, t1); //ar +  cr - dr
    Er = SUB(ar, t1); //ar - (cr - dr)
    
    Br = ADD(br, t0); // br + cr + dr
    Fr = SUB(t0, br); //-br + cr + dr
    
    c2 = LOAD(w2[0]);
    s2 = LOAD(w2[V]);

 
    //    REAL_TWISTED(Ar,Br,c2,s2,a4[0],a5[0]);

    real_twisted(Ar, Br, c2, s2, a4, a5);

    //c3 + i s3 = (c + i s) (c2 + i s2)
    /* t0 = MUL(c, c2); */
    /* t1 = MUL(s, s2); */
    /* t2 = MUL(c, s2); */
    /* t3 = MUL(s, c2); */

    /* c3 = SUB(t0, t1); */
    /* s3 = ADD(t2, t3); */

    t0 = mul0(c, c2);
    t1 = mul0(s, s2);
    t2 = mul0(c, s2);
    t3 = mul0(s, c2);

    c3 = SUB(t0, t1);
    s3 = ADD(t2, t3);

    red1(&c3);
    red1(&s3);

 
    //REAL_TWISTED(Er,Fr,c3,s3,a6[0],a7[0]);

    real_twisted(Er, Fr, c3, s3, a6, a7);
    
    a  += V;
    a1 += V;
    a2 += V;
    a3 += V;
    a4 += V;
    a5 += V;
    a6 += V;
    a7 += V;

    w  += V2;
    w2 += V2;
    
  }  while (k -= 1);
}


