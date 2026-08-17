#include "ntt.h"
#include "pre8_V4.h"

#include "i.h"


/* 
   FFT16V - 16*V complex numbers

   
   V <= 16

   a: 16*V complex = 32*V values  interleaved V-packed format  - V real parts, V imag parts, etc.

*/

void us16_vert(integer *a)
{

  //n = N/8 = 2

   const integer_packed p_64 = BROADCAST(PRIME_64);
   
  integer_packed x0, x1, x2, x3, x4, x5, x6, x7;
  integer_packed y0, y1, y2, y3, y4, y5, y6, y7;
  integer_packed xp04, xp26, xp15, xp37, yp04, yp26, yp15, yp37;
  //  integer_packed t02, s02, t13, s13;
  integer_packed c, s, Ar, Br, Ai, Bi, AipBr, AimBr, ArpBi, ArmBi;
  integer_packed Er, Ei, Fr, Fi;
  integer_packed ar,br,cr,dr,ai,bi,ci,di;
  integer_packed t1, t2, t4, t5;
  integer_packed s1, s2, s4, s5;
  integer_packed X0, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12, X13, X14, X15;
  integer_packed Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9, Y10, Y11, Y12, Y13, Y14, Y15;
  // integer_packed t0r, t0i, t1r, mt1i; //needed in macro UNTWISTED_REG

  /* z0 z1  z2  z3  z4  z5  z6  z7  z8  z9  z10  z11  z12  z13  z14  z15 */ 
  /*  0 2V  4V  6V  8V 10V 12V 14V 16V 18V  20V  22V  24V  26V  28V  30V */
  /*  0 V2 2V2 3V2 4V2 5V2 6V2 7V2 8V2 9V2 10V2 11V2 12V2 13V2 14V2 15V2 */


  //C4
  // U4_REG_MEM(a[0], a[V], a[2*V], a[3*V], a[4*V], a[5*V], a[6*V], a[7*V], X0, Y0, X1, Y1, X2, Y2, X3, Y3);

  u4_reg_to_mem(a+0*V, a+1*V, a+2*V, a+3*V, a+4*V, a+5*V, a+6*V, a+7*V, &X0, &Y0, &X1, &Y1, &X2, &Y2, &X3, &Y3);

  //C2 - 6 times

  /* U2_REG_MEM(a[8*V], a[9*V], a[10*V], a[11*V], X4, Y4, X5, Y5); */
  
  /* U2_REG_MEM(a[12*V], a[13*V], a[14*V], a[15*V], X6, Y6, X7, Y7); */

  /* U2_REG_MEM(a[16*V], a[17*V], a[18*V], a[19*V], X8, Y8, X9, Y9); */
  
  /* U2_REG_MEM(a[20*V], a[21*V], a[22*V], a[23*V], X10, Y10, X11, Y11); */
   
  /* U2_REG_MEM(a[24*V], a[25*V], a[26*V], a[27*V], X12, Y12, X13, Y13); */
    
  /* U2_REG_MEM(a[28*V], a[29*V], a[30*V], a[31*V], X14, Y14, X15, Y15); */

  
  u2_reg_to_mem(a+8*V, a+9*V, a+10*V, a+11*V, &X4, &Y4, &X5, &Y5);
  
  u2_reg_to_mem(a+12*V, a+13*V, a+14*V, a+15*V, &X6, &Y6, &X7, &Y7);

  u2_reg_to_mem(a+16*V, a+17*V, a+18*V, a+19*V, &X8, &Y8, &X9, &Y9);
  
  u2_reg_to_mem(a+20*V, a+21*V, a+22*V, a+23*V, &X10, &Y10, &X11, &Y11);
   
  u2_reg_to_mem(a+24*V, a+25*V, a+26*V, a+27*V, &X12, &Y12, &X13, &Y13);
    
  u2_reg_to_mem(a+28*V, a+29*V, a+30*V, a+31*V, &X14, &Y14, &X15, &Y15);

  

  //w = exp(2 pi i / 16)
  c = BROADCAST(1556715293UL);
  s = BROADCAST(978592373UL);
  //  c = BROADCAST(c_16[0]);
  // s = BROADCAST(s_16[0]);

  //w = exp(6 pi i / 16) = (s_16, c_16)

  //  UNTWISTED_REG(Er, Ei, Fr, Fi, s, c, X13, Y13, X11, Y11) w;

  untwisted_reg(&Er, &Ei, &Fr, &Fi, s, c, X13, Y13, X11, Y11);

 
  //  UNTWISTED_REG(Ar, Ai, Br, Bi, c, s, X9, Y9, X15, Y15);

  untwisted_reg(&Ar, &Ai, &Br, &Bi, c, s, X9, Y9, X15, Y15);

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
  

  //w0=1/sqrt(2)*(1+i), [ c=1, s=1]  then multiplication by sq2
  AimBr = ADD(X7, Y7);
  ArpBi = SUB(X7, Y7);

  ArmBi = ADD(Y5, X5);
  AipBr = SUB(Y5, X5);
  
  /* AipBr = MUL(q2, AipBr); */
  /* AimBr = MUL(q2, AimBr); */
    
  /* ArpBi = MUL(q2, ArpBi); */
  /* ArmBi = MUL(q2, ArmBi); */


  AipBr = ADD(AipBr, p_64);
  AipBr =  _mm256_slli_epi64(AipBr, SQ2_POW);
  red1(&AipBr);
  
  AimBr = ADD(AimBr, p_64);
  AimBr =  _mm256_slli_epi64(AimBr, SQ2_POW);
  red1(&AimBr);
  
  ArpBi = ADD(ArpBi, p_64);
  ArpBi =  _mm256_slli_epi64(ArpBi, SQ2_POW);
  red1(&ArpBi);
  
  ArmBi = ADD(ArmBi, p_64);
  ArmBi =  _mm256_slli_epi64(ArmBi, SQ2_POW);
  red1(&ArmBi);
  

  Ar = ADD(ArpBi, ArmBi);
  Bi = SUB(ArpBi, ArmBi);

  Ai = ADD(AipBr, AimBr);
  Br = SUB(AipBr, AimBr);

  xp04 = ADD(X1, Ar);
  xp26 = SUB(X1, Ar);

  xp15 = ADD(X3, Br);
  xp37 = SUB(X3, Br);

  x0 = ADD(xp04, ar);
  x4 = SUB(xp04, ar);

  x2 = ADD(xp26, br);
  x6 = SUB(xp26, br);

  x1 = ADD(xp15, cr);
  x5 = SUB(xp15, cr);

  x3 = ADD(xp37, dr);
  x7 = SUB(xp37, dr);

  yp04 = ADD(Y1, Ai);
  yp26 = SUB(Y1, Ai);
  
  yp15 = ADD(Y3, Bi);
  yp37 = SUB(Y3, Bi);
  
  y0 = ADD(yp04, ai);
  y4 = SUB(yp04, ai);

  y2 = ADD(yp26, bi);
  y6 = SUB(yp26, bi);

  y1 = ADD(yp15, ci);
  y5 = SUB(yp15, ci);

  y3 = ADD(yp37, di);
  y7 = SUB(yp37, di);

  STORE(a[1*V],  x0);
  STORE(a[3*V],  x1);
  STORE(a[9*V],  x2);
  STORE(a[11*V], x3);
  STORE(a[17*V], x4);
  STORE(a[19*V], x5);
  STORE(a[25*V], x6);
  STORE(a[27*V], x7);
  
  STORE(a[5*V],  y0);
  STORE(a[7*V],  y1);
  STORE(a[13*V], y2);
  STORE(a[15*V], y3);
  STORE(a[21*V], y4);
  STORE(a[23*V], y5);
  STORE(a[29*V], y6);
  STORE(a[31*V], y7);

  
  //finished second srl3 pass

  Ai = ADD(Y8, Y14);
  Br = SUB(Y8, Y14);

  Ar = ADD(X14, X8);
  Bi = SUB(X14, X8);
  
  Ei = ADD(Y12, Y10);
  Fr = SUB(Y12, Y10);

  Er = ADD(X10, X12);
  Fi = SUB(X10, X12);


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

  

  Ar = ADD(X6, X4);
  Bi = SUB(X6, X4);

  Ai = ADD(Y4, Y6);
  Br = SUB(Y4, Y6);


  xp04 = ADD(X0, Ar);
  xp26 = SUB(X0, Ar);

  xp15 = ADD(X2, Br);
  xp37 = SUB(X2, Br);

  x0 = ADD(xp04, ar);
  x4 = SUB(xp04, ar);

  x2 = ADD(xp26, br);
  x6 = SUB(xp26, br);

  x1 = ADD(xp15, cr);
  x5 = SUB(xp15, cr);

  x3 = ADD(xp37, dr);
  x7 = SUB(xp37, dr);

  yp04 = ADD(Y0, Ai);
  yp26 = SUB(Y0, Ai);
  
  yp15 = ADD(Y2, Bi);
  yp37 = SUB(Y2, Bi);
  
  y0 = ADD(yp04, ai);
  y4 = SUB(yp04, ai);

  y2 = ADD(yp26, bi);
  y6 = SUB(yp26, bi);

  y1 = ADD(yp15, ci);
  y5 = SUB(yp15, ci);

  y3 = ADD(yp37, di);
  y7 = SUB(yp37, di);

  STORE(a[0*V],   x0);
  STORE(a[2*V],   x1);
  STORE(a[8*V],   x2);
  STORE(a[10*V],  x3);
  STORE(a[16*V],  x4);
  STORE(a[18*V],  x5);
  STORE(a[24*V],  x6);
  STORE(a[26*V],  x7);
  
  STORE(a[4*V],  y0);
  STORE(a[6*V],  y1);
  STORE(a[12*V], y2);
  STORE(a[14*V], y3);
  STORE(a[20*V], y4);
  STORE(a[22*V], y5);
  STORE(a[28*V], y6);
  STORE(a[30*V], y7);


}
