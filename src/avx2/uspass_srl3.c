#include "ntt.h"

#include "pre8_V4.h"

#include "i.h"


void uspass_srl3(integer *a, const integer *w, const integer *w2, size_t n)
{

  size_t k;
    
  integer *a1;
  integer *a2;
  integer *a3;
  integer *a4;
  integer *a5;
  integer *a6;
  integer *a7;

  
  integer *b;
  integer *b1;
  integer *b2;
  integer *b3;
  integer *b4;
  integer *b5;
  integer *b6;
  integer *b7;


  const integer_packed p_64 = BROADCAST(PRIME_64);
  
  integer_packed x0, x1, x2, x3, x4, x5, x6, x7;
  integer_packed y0, y1, y2, y3, y4, y5, y6, y7;

  integer_packed A0r, B0r, A0i, B0i;
  integer_packed Ar, Br, Ai, Bi;
  integer_packed Er, Ei, Fr, Fi;

  integer_packed ar,br,cr,dr,ai,bi,ci,di;
 
  integer_packed t0, t1, t2, t3, t4, t5;
  integer_packed s0, s1, s2, s3, s4, s5;
  integer_packed c, s, cc2, ss2, cc3, ss3;
   
  
  
  a1 = a  + n;
  a2 = a1 + n;
  a3 = a2 + n;
  a4 = a3 + n;
  a5 = a4 + n;
  a6 = a5 + n;
  a7 = a6 + n;

  b  = a  + (n<<3);
  b1 = b  + n;
  b2 = b1 + n;
  b3 = b2 + n;
  b4 = b3 + n;
  b5 = b4 + n;
  b6 = b5 + n;
  b7 = b6 + n;
  
  
  k = n >> VPOW;  //n/V


  do {

    cc2 = LOAD(w2[0]);
    ss2 = LOAD(w2[V]);

    //    UNTWISTED(Ar, Ai, Br, Bi, cc2, ss2, b[0],b1[0],b6[0],b7[0]);

    untwisted(&Ar, &Ai, &Br, &Bi, cc2, ss2, b, b1, b6, b7);

    c = LOAD(w[0]);
    s = LOAD(w[V]);
    
    t0 = mul0(c, cc2);
    t1 = mul0(s, ss2);
    t2 = mul0(c, ss2);
    t3 = mul0(s, cc2);

    cc3 = SUB(t0, t1);
    ss3 = ADD(t2, t3);

    red1(&cc3);
    red1(&ss3);

    //    UNTWISTED(Er, Ei, Fr, Fi, cc3, ss3, b4[0],b5[0],b2[0],b3[0]);

    untwisted(&Er, &Ei, &Fr, &Fi, cc3, ss3, b4, b5, b2, b3);

    
    ar = ADD(Ar, Er); //A + E = a             X 0-4
    t1 = SUB(Ar, Er); //A - E = c_bar

    t2 = ADD(Br, Fr); //B + F = d_bar
    br = SUB(Br, Fr); //B - F = b             X 2-6

    t4 = ADD(t1, t2); //c_bar + d_bar = c     X 1-5
    t5 = SUB(t2, t1); //d_bar - c_bar = d     X 3-7

    ai = ADD(Ai, Ei); //A + E = a             Y 0-4
    s1 = SUB(Ai, Ei); //A - E = c_bar

    s2 = ADD(Bi, Fi); //B + F = d_bar
    bi = SUB(Bi, Fi); //B - F = b             Y 2-6

    s4 = ADD(s1, s2); //c_bar + d_bar = c     Y 1-5
    s5 = SUB(s2, s1); //d_bar - c_bar = d     Y 3-7
    

    /* q2 = BROADCAST(sq2[0]); */

    /* cr = MUL(q2, t4); //X 1-5 */
    /* dr = MUL(q2, t5); //X 3-7 */
    /* ci = MUL(q2, s4); //Y 1-5 */
    /* di = MUL(q2, s5); //Y 3-7 */

    cr = ADD(t4, p_64);
    cr = _mm256_slli_epi64(cr, SQ2_POW);
    red1(&cr);
    
    dr = ADD(t5, p_64);
    dr = _mm256_slli_epi64(dr, SQ2_POW);
    red1(&dr);
    
    ci = ADD(s4, p_64);
    ci = _mm256_slli_epi64(ci, SQ2_POW);
    red1(&ci);
    
    di = ADD(s5, p_64);
    di = _mm256_slli_epi64(di, SQ2_POW);
    red1(&di);

    
    //    UNTWISTED(A0r, A0i, B0r, B0i, c, s, a4[0],a5[0],a6[0],a7[0]); //0+4-2(2+6) integer and imag and 1+5-(3+7) real and imag

    untwisted(&A0r, &A0i, &B0r, &B0i, c, s, a4, a5, a6, a7);

    t0 = LOAD(a[0]); // X 0+4+2+6
    s0 = LOAD(a2[0]); // Y 0+4+2+6
    
    t1 = LOAD(a1[0]); // X 1+5+3+7
    s1 = LOAD(a3[0]); // Y 1+5+3+7


    
    t2 = ADD(t0, A0r); //X 0+4
    t3 = SUB(t0, A0r); //X 2+6

    s2 = ADD(s0, A0i); //Y 0+4
    s3 = SUB(s0, A0i); //Y 2+6

    x0 = ADD(t2, ar);
    x4 = SUB(t2, ar);

    y0 = ADD(s2, ai);
    y4 = SUB(s2, ai);

    STORE(a[0], x0);
    STORE(b[0], y0);
    
    STORE(a4[0], x4);
    STORE(b4[0], y4);

    x2 = ADD(t3, br);
    x6 = SUB(t3, br);

    y2 = ADD(s3, bi);
    y6 = SUB(s3, bi);

    STORE(a2[0], x2);
    STORE(b2[0], y2);
    
    STORE(a6[0], x6);
    STORE(b6[0], y6);

    
    t2 = ADD(t1, B0r); //X 1+5
    t3 = SUB(t1, B0r); //X 3+7

    s2 = ADD(s1, B0i); //Y 1+5
    s3 = SUB(s1, B0i); //Y 3+7
       
    x1 = ADD(t2, cr);
    x5 = SUB(t2, cr);

    y1 = ADD(s2, ci);
    y5 = SUB(s2, ci);

    STORE(a1[0], x1);
    STORE(b1[0], y1);
    
    STORE(a5[0], x5);
    STORE(b5[0], y5);

   
    x3 = ADD(t3, dr);
    x7 = SUB(t3, dr);

    y3 = ADD(s3, di);
    y7 = SUB(s3, di);

    STORE(a3[0], x3);
    STORE(b3[0], y3);
    
    STORE(a7[0], x7);
    STORE(b7[0], y7);

    
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
