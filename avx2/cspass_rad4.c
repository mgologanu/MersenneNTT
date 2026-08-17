#include "ntt.h"
#include "pre8_V4.h"

#include "i.h"

//#include <stdio.h>

/* 
   Radix 4 pass for x^8n-1 - 8n complex numbers 

   a: 8n complex = 16n values = split real and imag arrays

   w = w0, w
   w0: 1/2th roots of order 4n
   w:  1/4th roots of order 8n



*/
void cspass_rad4(integer *a,  const integer *w,  size_t n)
{

  
  /* const integer_packed p = BROADCAST(PRIME); */

  /* const integer_packed p_64 = BROADCAST(PRIME_64); */

  
  // integer_packed _mul1_t1, _mul1_t2, _mul1_t3, retval;
  //integer_packed _red1_t1, _red1_t2, _red1_t3;
  
  // integer_packed _t0, _t1, _t2, _AipBr, _AimBr, _ArpBi, _ArmBi;
  
  integer_packed x0, x1, x2, x3, y0, y1, y2, y3;
  integer_packed c, s, t02, t13, s02, s13, Ar, Br, Ai, Bi, tt0, tt1, tt2, tt3;
  integer_packed Xm01, Ym01;
  
  size_t k, n2;
  
  
  integer *a1;
  integer *a2;
  integer *a3;

  integer *b0;
  integer *b1;
  integer *b2;
  integer *b3;
  

  n2 = n << 1; //2n
  
  a1 = a  + n2;
  a2 = a1 + n2;
  a3 = a2 + n2;

  b0 = a3 + n2;
  b1 = b0 + n2;
  b2 = b1 + n2;
  b3 = b2 + n2;
  
  

  k = n >> (VPOW-1);  //2n/V
  //  printf("k = %d, n2 = %d\n", k, n2);
  do {
    
    x0 = LOAD(a[0]);
    x1 = LOAD(a1[0]);
    x2 = LOAD(a2[0]);
    x3 = LOAD(a3[0]);
    
    y0 = LOAD(b0[0]);
    y1 = LOAD(b1[0]);
    y2 = LOAD(b2[0]);
    y3 = LOAD(b3[0]);

    t02 = ADD(x0, x2);
    Ar  = SUB(x0, x2);
    
    t13 = ADD(x1, x3);
    Br  = SUB(x1, x3);
    
    s02 = ADD(y0, y2);
    Ai  = SUB(y0, y2);
    
    s13 = ADD(y1, y3);
    Bi  = SUB(y1, y3);

    tt0 = ADD(t02, t13); //X 0+2+1+3
    tt1 = ADD(s02, s13); //Y 0+2+1+3

    STORE(a[0], tt0);
    STORE(a1[0], tt1);
    
    Xm01 = SUB(t02, t13); //X 0+2-1-3
    Ym01=  SUB(s02, s13); //Y 0+2-1-3

    //(Xm01 + i * Ym01) * (C0 + i * S0)
    c = LOAD(w[0]);//C0
    s = LOAD(w[V]);//S0

    tt0 = mul1(Xm01, c);
    tt1 = mul1(Ym01, s);
    tt2 = mul1(Ym01, c);
    tt3 = mul1(Xm01, s);
    
    tt0 = SUB(tt0, tt1);
    tt2 = ADD(tt2, tt3);

    STORE(a2[0], tt0);
    STORE(a3[0], tt2);
    
    c = LOAD(w[2*V]);//C
    s = LOAD(w[3*V]);//S
    
    //    TWISTED(Ar, Ai, Br, Bi, c, s, b0[0], b1[0], b2[0], b3[0]);

    twisted(Ar, Ai, Br, Bi, c, s, b0, b1, b2, b3);
    a  += V;
    a1 += V;
    a2 += V;
    a3 += V;

    b0 += V;
    b1 += V;
    b2 += V;
    b3 += V;

    w  += (V<<3); //4V - 2(cos, sin) = w0, w, 8V if using rad8 roots
    
  }  while (k -= 1);
}



