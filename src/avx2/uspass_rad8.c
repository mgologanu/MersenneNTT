#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void uspass_rad8(integer *a,  const integer *w, size_t n)
{

  size_t k;
  const integer_packed p_64 = BROADCAST(PRIME_64);

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
  
  integer_packed x0, x1, x2, x3, x4, x5, x6, x7;
  integer_packed y0, y1, y2, y3, y4, y5, y6, y7;

  integer_packed A0r, B0r, A0i, B0i;
  integer_packed Ar, Br, Ai, Bi;
  integer_packed Er, Ei, Fr, Fi;
  integer_packed X0, X1, Y0, Y1, Xm01, Ym01; 

  integer_packed ar,br,cr,dr,ai,bi,ci,di;
  
  integer_packed t0, t1, t2, t3, t4, t5;
  integer_packed s0, s1, s2, s3, s4, s5;
  integer_packed c, s;
 
  
  
  k = n >> VPOW;  //n/V


  do {
    c = LOAD(w[4*V]);
    s = LOAD(w[5*V]);

    //    UNTWISTED(Ar, Ai, Br, Bi, c, s, b[0],b1[0],b6[0],b7[0]);
    
    untwisted(&Ar, &Ai, &Br, &Bi, c, s, b, b1, b6, b7);
 
    c = LOAD(w[6*V]);
    s = LOAD(w[7*V]);


    //    UNTWISTED(Er, Ei, Fr, Fi, c, s, b4[0],b5[0],b2[0],b3[0]);

    untwisted(&Er, &Ei, &Fr, &Fi, c, s, b4, b5, b2, b3);

    
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

    
    c = LOAD(w[2*V]);//C
    s = LOAD(w[3*V]);//S
   
    //UNTWISTED(A0r, A0i, B0r, B0i, c, s, a4[0],a5[0],a6[0],a7[0]); //0+4-2(2+6) integer and imag and 1+5-(3+7) real and imag
    untwisted(&A0r, &A0i, &B0r, &B0i, c, s, a4, a5, a6, a7);


    c = LOAD(w[0]);//C0
    s = LOAD(w[V]);//S0

    X1 = LOAD(a2[0]);
    Y1 = LOAD(a3[0]);

    //X1 + i* Y1 = (Xm01 + i * Ym01) * (C0 + i * S0)  --> Xm01 + i * Ym01 = (X1 + i * Y1)*(C0 - i * S0) = X1 * C0 + Y 1 * S0 + i (Y1 C0 - X1 S0)
    
    t0 = mul1(X1, c);
    t1 = mul1(Y1, s);
    t2 = mul1(Y1, c);
    t3 = mul1(X1, s);

    Xm01 = ADD(t0, t1);  //X 0+4+2+6-(1+5+3+7)
    Ym01 = SUB(t2, t3);  //Y 0+4+2+6-(1+5+3+7)

    X0 = LOAD(a[0]);     //X 0+4+2+6+(1+5+3+7)
    Y0 = LOAD(a1[0]);     //Y 0+4+2+6+(1+5+3+7)


    t0 = ADD(X0, Xm01); //X 0+4+2+6
    s0 = ADD(Y0, Ym01); //Y 0+4+2+6

    
    t1 = SUB(X0, Xm01); //X 1+5+3+7
    s1 = SUB(Y0, Ym01); //Y 1+5+3+7

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


    w  += (V<<3); //8V - 4 pairs (cos, sin) = w0, w, w2, w3
  
    
  }  while (k -= 1);

    
}
