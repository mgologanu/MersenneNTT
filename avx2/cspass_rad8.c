#include "ntt.h"
#include "pre8_V4.h"

#include "i.h"

/* 
   Radix 8  pass for x^8n-1 - 8n complex numbers 

   a: 8n complex = 16n values in split real and imag arrays, but contiguous

   w: 4n complex numbers or 8n values: 1/2th roots of 2n, 1/4th roots
   of order 4n, 1/8th roots of order 8n, 1/8th roots of order 8n at power 3
      
   w = w0, w, w2, w3
*/
void cspass_rad8(integer *a,  const integer *w, size_t n)
{
  //  const integer_packed  p = BROADCAST(PRIME);

    const integer_packed p_64 = BROADCAST(PRIME_64);

  
  // integer_packed _mul1_t1, _mul1_t2, _mul1_t3, retval;
  // integer_packed _red1_t1, _red1_t2, _red1_t3;
  
  //  integer_packed _t0, _t1, _t2,_AipBr, _AimBr, _ArpBi, _ArmBi;
  
    
  integer_packed x0, x1, x2, x3, x4, x5, x6, x7;
  integer_packed y0, y1, y2, y3, y4, y5, y6, y7;
  integer_packed xp04, xp26, xp15, xp37, yp04, yp26, yp15, yp37; 
  integer_packed c, s, Ar, Br, Ai, Bi, t0, t1, t2, t3;
  integer_packed Er, Ei, Fr, Fi;
  integer_packed ar,br,cr,dr,ai,bi,ci,di;
  integer_packed X0,X1,Y0,Y1, Xm01, Ym01;

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

    X0 = ADD(xp04, xp26);
    Ar = SUB(xp04, xp26);
 

    X1 = ADD(xp15, xp37);
    Br = SUB(xp15, xp37);

    t0 = ADD(X0, X1); //X 0+4+2+6+1+5+3+7 
    STORE(a[0], t0);  //X0

    Xm01 = SUB(X0, X1);  //X 0+4+2+6-(1+5+3+7 )

    yp04 = ADD(y0, y4);
    ai   = SUB(y0, y4);
    
    yp26 = ADD(y2, y6);
    bi   = SUB(y2, y6);

    yp15 = ADD(y1, y5);
    ci   = SUB(y1, y5);
    
    yp37 = ADD(y3, y7);
    di   = SUB(y3, y7);

    Y0 = ADD(yp04, yp26);
    Ai = SUB(yp04, yp26);


    Y1 = ADD(yp15, yp37); //Y 0+4+2+6+1+5+3+7 
    Bi = SUB(yp15, yp37);

    t0 = ADD(Y0, Y1);
    STORE(a1[0], t0);  //Y0

    Ym01 = SUB(Y0, Y1); //Y 0+4+2+6-(1+5+3+7 )

    //(Xm01 + i * Ym01) * (C0 + i * S0)
    c = LOAD(w[0]);//C0
    s = LOAD(w[V]);//S0
    
    t0 = mul1(Xm01, c);
    t1 = mul1(Ym01, s);
    t2 = mul1(Ym01, c);
    t3 = mul1(Xm01, s);

    t0 = SUB(t0, t1);
    t2 = ADD(t2, t3);

    STORE(a2[0], t0);//X1
    STORE(a3[0], t2);//Y1
    
    c = LOAD(w[2*V]);//C
    s = LOAD(w[3*V]);//S

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
    

    t0 = ADD(cr, dr);
    t1 = SUB(cr, dr);

    Ar = ADD(ar, t1);
    Er = SUB(ar, t1);
    
    Br = ADD(br, t0);
    Fr = SUB(t0, br);
    
    t0 = ADD(ci, di);
    t1 = SUB(ci, di);

    Ai = ADD(ai, t1);
    Ei = SUB(ai, t1);

    
    Bi = ADD(bi, t0);
    Fi = SUB(t0, bi);

    c = LOAD(w[4*V]);//C2
    s = LOAD(w[5*V]);//S2

    //    TWISTED(Ar,Ai,Br,Bi,c,s,b[0],b1[0],b6[0],b7[0]);

    twisted(Ar, Ai, Br, Bi, c, s, b, b1, b6, b7);

    c = LOAD(w[6*V]);//C3
    s = LOAD(w[7*V]);//S3

    //    TWISTED(Er,Ei,Fr,Fi,c,s,b4[0],b5[0],b2[0],b3[0]);

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

    w  += (V8);  //8V - 4 x (cos, sin) = w0, w, w2, w3
    
  }  while (k -= 1);
}

