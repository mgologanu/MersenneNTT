#include "ntt.h"
#include "pre8_V4.h"

#include "i.h"

/* 
   FFT8V - 8*V complex numbers

   
   V <= 8

   a: 8*V complex = 16*V values  interleaved V-packed format  - V real parts, V imag parts, etc.

*/
void c8_vert(integer *a)
{

  //n = N/8 = 1
  const integer_packed  p = BROADCAST(PRIME);
  
  const integer_packed p_64 = BROADCAST(PRIME_64);
  
  integer_packed _t0, _t1, _t2, _t3, retval;

  integer_packed _red1_t1, _red1_t2, _red1_t3;
  
  integer_packed x0, x1, x2, x3, x4, x5, x6, x7;
  integer_packed y0, y1, y2, y3, y4, y5, y6, y7;
  integer_packed xp04, xp26, xp15, xp37, yp04, yp26, yp15, yp37; 
  integer_packed Ar, Br, Ai, Bi, AipBr, AimBr, ArpBi, ArmBi, t0, t1;
  integer_packed Er, Ei, Fr, Fi;
  integer_packed ar,br,cr,dr,ai,bi,ci,di;
  integer_packed X0,X1,Y0,Y1, Xm01, Ym01;
 
  /* 
        ---C8 -----  ---C8 -----  ---C8 ----- ---C8 ----- 
        x0 y0 x1 y1  x2 y2 x3 y3  x4 y4 x5 y5  x6 y6 x7 y7 
    (V*) 0  1  2  3   4  5  6  7   8  9 10 11  12 13 14 15
     
    transpose(x0 x2 x4 x6) -> t0 t1 t2 t3 
    
    transpose(x1 x3 x5 x7) -> t4 t5 t6 t7
    
    save in place in increasing order of ti
    
    Repeat for y
      
  */
  

  /* x0 = LOAD(a[0]); */
  /* x1 = LOAD(a[4*V]); */
  /* x2 = LOAD(a[8*V]); */
  /* x3 = LOAD(a[12*V]); */

  /* TRANSPOSE4(x0,x1,x2,x3); */

  TRANSPOSE4_FROM_MEM(x0,x1,x2,x3,a,0,V4);
  
  /* x4 = LOAD(a[2*V]); */
  /* x5 = LOAD(a[6*V]); */
  /* x6 = LOAD(a[10*V]); */
  /* x7 = LOAD(a[14*V]); */
  
  /* TRANSPOSE4(x4,x5,x6,x7); */

  TRANSPOSE4_FROM_MEM(x4,x5,x6,x7,a, V2, V4);


  /* y0 = LOAD(a[1*V]); */
  /* y1 = LOAD(a[5*V]); */
  /* y2 = LOAD(a[9*V]); */
  /* y3 = LOAD(a[13*V]); */

  /* TRANSPOSE4(y0,y1,y2,y3); */

  TRANSPOSE4_FROM_MEM(y0,y1,y2,y3, a, V, V4);

  /* y4 = LOAD(a[3*V]); */
  /* y5 = LOAD(a[7*V]); */
  /* y6 = LOAD(a[11*V]); */
  /* y7 = LOAD(a[15*V]); */
  
  /* TRANSPOSE4(y4,y5,y6,y7); */

  TRANSPOSE4_FROM_MEM(y4,y5,y6,y7, a, 3*V, V4);
  
    
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
  
  t0 = ADD(X0, X1);
  STORE(a[0], t0);  
  
  Xm01 = SUB(X0, X1);
  
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
  
  
  Y1 = ADD(yp15, yp37);
  Bi = SUB(yp15, yp37);
  
  t0 = ADD(Y0, Y1);
  STORE(a[V], t0);
  
  Ym01 = SUB(Y0, Y1);

  
  // c = 1, s = 0
  
  STORE(a[2*V], Xm01);  
  STORE(a[3*V], Ym01);

  // c = 1, s = 0

  AipBr = ADD(Ai, Br);
  AimBr = SUB(Ai, Br);
  
  ArpBi = ADD(Ar, Bi);
  ArmBi = SUB(Ar, Bi);
  
  STORE(a[4*V], ArmBi); 
  STORE(a[5*V], AipBr);

  STORE(a[6*V], ArpBi);  
  STORE(a[7*V], AimBr);
 
  /* q2 = BROADCAST(sq2[0]); */
  /* cr = MUL(q2, cr); */
  /* dr = MUL(q2, dr); */
  /* ci = MUL(q2, ci); */
  /* di = MUL(q2, di); */

  //Multiplication by 2^15, as sqrt(i) = 2^15*(1+i)

  cr = RED1(cr);
  cr =  _mm256_slli_epi64(cr, SQ2_POW);
  cr = RED1(cr);

  dr = RED1(dr);
  dr =  _mm256_slli_epi64(dr, SQ2_POW);
  dr = RED1(dr);

  ci = RED1(ci);
  ci =  _mm256_slli_epi64(ci, SQ2_POW);
  ci = RED1(ci);

  di = RED1(di);
  di =  _mm256_slli_epi64(di, SQ2_POW);
  di = RED1(di);
  
    
  
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

  
  //c = 1, s=0
   
  AipBr = ADD(Ai, Br);
  AimBr = SUB(Ai, Br);
  
  ArpBi = ADD(Ar, Bi);
  ArmBi = SUB(Ar, Bi);
  
  STORE(a[8*V], ArmBi);  
  STORE(a[9*V], AipBr);

  STORE(a[14*V], ArpBi);  
  STORE(a[15*V], AimBr);

  //c = 1, s = 0

  AipBr = ADD(Ei, Fr);
  AimBr = SUB(Ei, Fr);
  
  ArpBi = ADD(Er, Fi);
  ArmBi = SUB(Er, Fi);
  
  STORE(a[12*V], ArmBi);    
  STORE(a[13*V], AipBr);

  STORE(a[10*V], ArpBi);     
  STORE(a[11*V], AimBr); 
  
}


