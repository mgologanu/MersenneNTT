#include "ntt.h"
#include "pre8_V4.h"

#include "i.h"

/* 
   FFT8V - 8*V complex numbers

   
   V <= 8

   a: 8*V complex = 16*V values in a partially split format

*/
void cs8_vert(integer *a)
{

  //n = N/8 = 1
  /* const integer_packed  p = BROADCAST(PRIME); */
  
   const integer_packed p_64 = BROADCAST(PRIME_64); 
  
  /* integer_packed _t0, _t1, _t2, _t3; */

  /* integer_packed _red1_t1, _red1_t2, _red1_t3, retval; */
  
  __m256i x0, x1, x2, x3, x4, x5, x6, x7;
  __m256i y0, y1, y2, y3, y4, y5, y6, y7;
  __m256i xp04, xp26, xp15, xp37, yp04, yp26, yp15, yp37; 
  __m256i Ar, Br, Ai, Bi, AipBr, AimBr, ArpBi, ArmBi, t0, t1;
  __m256i Er, Ei, Fr, Fi;
  __m256i ar,br,cr,dr,ai,bi,ci,di;
  __m256i X0,X1,Y0,Y1, Xm01, Ym01;

  /* integer_packed tt1, tt2, tt3; */

  /* integer va[4]; */
  
  /*

    This is a partially split format from the output of cspass_rad8

        ---C8 -----  ---C8 -----  ---C8 -----   ---C8 -----   
        x0 x1 y0 y1  x2 x3 y2 y3  x4 x5 y4 y5   x6 x7 y6 y7
    (V*) 0  1  2  3   4  5  6  7   8  9 10 11   12 13 14 15 
    
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

  // TRANSPOSE4_FROM_MEM(x0,x1,x2,x3,a,0,V4);

  /* x4 = LOAD(a[1*V]); */
  /* x5 = LOAD(a[5*V]); */
  /* x6 = LOAD(a[9*V]); */
  /* x7 = LOAD(a[13*V]); */
  
  /* TRANSPOSE4(x4,x5,x6,x7); */

  //  TRANSPOSE4_FROM_MEM(x4,x5,x6,x7,a, V, V4);
  
  /* y0 = LOAD(a[2*V]); */
  /* y1 = LOAD(a[6*V]); */
  /* y2 = LOAD(a[10*V]); */
  /* y3 = LOAD(a[14*V]); */

  /* TRANSPOSE4(y0,y1,y2,y3); */
  
  //  TRANSPOSE4_FROM_MEM(y0,y1,y2,y3, a, V2, V4);
   
  /* y4 = LOAD(a[3*V]); */
  /* y5 = LOAD(a[7*V]); */
  /* y6 = LOAD(a[11*V]); */
  /* y7 = LOAD(a[15*V]); */

  /* TRANSPOSE4(y4,y5,y6,y7); */

  //  TRANSPOSE4_FROM_MEM(y4,y5,y6,y7, a, 3*V, V4);


  transpose4_from_mem(&x0, &x1, &x2, &x3, a, 0, V4);
  transpose4_from_mem(&x4, &x5, &x6, &x7, a, V, V4);
  transpose4_from_mem(&y0, &y1, &y2, &y3, a, V2, V4);
  transpose4_from_mem(&y4, &y5, &y6, &y7, a, 3*V, V4);
  

  /* Do C8 for V-packed numbers */
    
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

  /* q2 = BROADCAST(32768UL);  */
  /* cr = MUL1(cr, q2); */
  /* dr = MUL1(dr, q2); */
  /* ci = MUL1(ci, q2); */
  /* di = MUL1(di, q2); */

  
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


