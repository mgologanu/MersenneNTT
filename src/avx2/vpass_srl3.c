#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"


void vpass_srl3(integer *a, const integer *w, const integer *w2, size_t n)
{

  const integer_packed p_64 = BROADCAST(PRIME_64);
  
  integer_packed x0, x1, x2, x3, x4, x5, x6, x7;
  integer_packed c, s, c2, s2, c3, s3;
  integer_packed Ar, Br, A0r, B0r;
  integer_packed t0, t1, t2, t3, t4, t5;
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

    c2 = LOAD(w2[0]);
    s2 = LOAD(w2[V]);

    //REAL_UNTWISTED(Ar,Br,c2,s2,a4[0],a5[0]);
    
    real_untwisted(&Ar, &Br, c2, s2, a4, a5);
    
    c = LOAD(w[0]);
    s = LOAD(w[V]);
    
    t0 = mul0(c, c2);
    t1 = mul0(s, s2);
    t2 = mul0(c, s2);
    t3 = mul0(s, c2);

    c3 = SUB(t0, t1);
    s3 = ADD(t2, t3);
    
    red1(&c3);
    red1(&s3);

    //    INTEGER_UNTWISTED(Er,Fr,c3,s3,a6[0],a7[0]);
    real_untwisted(&Er, &Fr, c3, s3, a6, a7);
 
    

    ar = ADD(Ar, Er); //A + E = a             X 0-4
    t1 = SUB(Ar, Er); //A - E = c_bar

    t2 = ADD(Br, Fr); //B + F = d_bar
    br = SUB(Br, Fr); //B - F = b             X 2-6

    t4 = ADD(t1, t2); //c_bar + d_bar = c     X 1-5
    t5 = SUB(t2, t1); //d_bar - c_bar = d     X 3-7

   
    /* q2 = BROADCAST(sq2[0]); */

    /* cr = MUL(q2, t4); //X 1-5 */
    /* dr = MUL(q2, t5); //X 3-7 */

    cr = ADD(t4, p_64);
    cr = _mm256_slli_epi64(cr, SQ2_POW);
    red1(&cr);
    
    dr = ADD(t5, p_64);
    dr = _mm256_slli_epi64(dr, SQ2_POW);
    red1(&dr);
    

    // REAL_UNTWISTED(A0r,B0r,c,s,a2[0],a3[0]);
    real_untwisted(&A0r, &B0r, c, s, a2, a3);

     
    t0 = LOAD(a[0]); // X 0+4+2+6
    t1 = LOAD(a1[0]); // X 1+5+3+7

    t2 = ADD(t0, A0r); //X 0+4
    t3 = SUB(t0, A0r); //X 2+6

    x0 = ADD(t2, ar);
    x4 = SUB(t2, ar);
    STORE(a[0],  x0);
    STORE(a4[0], x4);
 

    x2 = ADD(t3, br);
    x6 = SUB(t3, br);
    STORE(a2[0], x2);    
    STORE(a6[0], x6);


    t2 = ADD(t1, B0r); //X 1+5
    t3 = SUB(t1, B0r); //X 3+7
    
    x1 = ADD(t2, cr);
    x5 = SUB(t2, cr);
    STORE(a1[0], x1);    
    STORE(a5[0], x5);
   
    x3 = ADD(t3, dr);
    x7 = SUB(t3, dr);
    STORE(a3[0], x3);
    STORE(a7[0], x7);

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
