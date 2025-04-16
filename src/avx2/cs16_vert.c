#include "ntt.h"
#include "pre8_V4.h"

#include "i.h"


/* 
   FFT16V - 16*V complex numbers

   
   V <= 16

   a: 16*V complex = 32*V values as partially split real and imag arrays, given by the output of rad8 (same as t16):


   x0 x1 x2 x3  y0 y1 y2 y3  x4 x5 x6 x7  y4 y5 y6 y7  x8 x9 x10 x11  y8 y9 y10 y11  x12 x13 x14 x15  y12 y13 y14 y15 
    0  1  2  3  4  5  6  7   8  9 10  11  12 13 14 15  16 17 18  19  20 21  22  23   24  25  26  27   28  29  30  31   


*/

void cs16_vert(integer *a)
{

  //n = N/8 = 2
  /* const integer_packed  p = BROADCAST(PRIME); */
  
  const integer_packed p_64 = BROADCAST(PRIME_64);

  /* integer_packed _mul1_t1, _mul1_t2, _mul1_t3, retval; */
  /* integer_packed _red1_t1, _red1_t2, _red1_t3; */
  
  /* integer_packed _t0, _t1, _t2, _t3, _AipBr, _AimBr, _ArpBi, _ArmBi; */

  /* integer_packed _Ar, _Br, _Ai, _Bi, _t02, _t13, _s02, _s13; */

  integer_packed AipBr, AimBr, ArpBi, ArmBi;

  integer_packed x0, x1, x2, x3, x4, x5, x6, x7;
  integer_packed y0, y1, y2, y3, y4, y5, y6, y7;
  integer_packed xp04, xp26, xp15, xp37, yp04, yp26, yp15, yp37;
  integer_packed c, s, Ar, Br, Ai, Bi;
  integer_packed Er, Ei, Fr, Fi;
  integer_packed ar,br,cr,dr,ai,bi,ci,di;
  integer_packed t0, t1;
  integer_packed X0, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11, X12, X13, X14, X15;
  integer_packed Y0, Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9, Y10, Y11, Y12, Y13, Y14, Y15;

  //First srl3 pass, all c,s = 1,0
  // x0 x2 x4 x6 x8 x10 x12 x14
  // same for y

  x0 = LOAD(a[0*V]);
  x1 = LOAD(a[2*V]);
  
  x2 = LOAD(a[8*V]);
  x3 = LOAD(a[10*V]);

  x4 = LOAD(a[16*V]);
  x5 = LOAD(a[18*V]);
  
  x6 = LOAD(a[24*V]);
  x7 = LOAD(a[26*V]);

  
  y0 = LOAD(a[4*V]);
  y1 = LOAD(a[6*V]);
  
  y2 = LOAD(a[12*V]);
  y3 = LOAD(a[14*V]);
  
  y4 = LOAD(a[20*V]);
  y5 = LOAD(a[22*V]);

  
  y6 = LOAD(a[28*V]);
  y7 = LOAD(a[30*V]);

  
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
   
  X2 = ADD(xp15, xp37);
  Br = SUB(xp15, xp37);
   
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
   
  Y2 = ADD(yp15, yp37);
  Bi = SUB(yp15, yp37);
   
  // c = 1, s = 0
  
  X4 = SUB(Ar, Bi);
  Y4 = ADD(Ai, Br);
  X6 = ADD(Ar, Bi);
  Y6 = SUB(Ai, Br);

  /* q2 = BROADCAST(32768UL);  */
  /* /\* q2 = BROADCAST(sq2[0]); *\/ */
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

  
  // c = 1, s=0

  X8  = SUB(Ar, Bi);
  Y8  = ADD(Ai, Br);
  X14 = ADD(Ar, Bi);
  Y14 = SUB(Ai, Br);
    
  //c = 1, s = 0
 
  X12 = SUB(Er, Fi);
  Y12 = ADD(Ei, Fr);
  X10 = ADD(Er, Fi);
  Y10 = SUB(Ei, Fr);

  //second srl3 pass

  //x1 x3 x5 x7 x9 x11 x13 x15
  
  x0 = LOAD(a[1*V]);
  x1 = LOAD(a[3*V]);
  
  x2 = LOAD(a[9*V]);
  x3 = LOAD(a[11*V]);

  x4 = LOAD(a[17*V]);
  x5 = LOAD(a[19*V]);

  x6 = LOAD(a[25*V]);
  x7 = LOAD(a[27*V]);
  
  y0 = LOAD(a[5*V]);
  y1 = LOAD(a[7*V]);
  
  y2 = LOAD(a[13*V]);
  y3 = LOAD(a[15*V]);
  
  y4 = LOAD(a[21*V]);
  y5 = LOAD(a[23*V]);

  
  y6 = LOAD(a[29*V]);
  y7 = LOAD(a[31*V]);


  xp04 = ADD(x0, x4);
  ar   = SUB(x0, x4);
    
  xp26 = ADD(x2, x6);
  br   = SUB(x2, x6);

  xp15 = ADD(x1, x5);
  cr   = SUB(x1, x5);
    
  xp37 = ADD(x3, x7);
  dr   = SUB(x3, x7);

  X1 = ADD(xp04, xp26);
  Ar = SUB(xp04, xp26);
 
  X3 = ADD(xp15, xp37);
  Br = SUB(xp15, xp37);
  
  yp04 = ADD(y0, y4);
  ai   = SUB(y0, y4);
    
  yp26 = ADD(y2, y6);
  bi   = SUB(y2, y6);

  yp15 = ADD(y1, y5);
  ci   = SUB(y1, y5);
    
  yp37 = ADD(y3, y7);
  di   = SUB(y3, y7);

  Y1 = ADD(yp04, yp26);
  Ai = SUB(yp04, yp26);
 
  Y3 = ADD(yp15, yp37);
  Bi = SUB(yp15, yp37);
  
  //w0=1/sqrt(2)*(1+i),  c=1, s=1 + multiplication by sq2
    
  AipBr = ADD(Ai, Br);
  AimBr = SUB(Ai, Br);
    
  ArpBi = ADD(Ar, Bi);
  ArmBi = SUB(Ar, Bi);

        
  /* AipBr = MUL1(AipBr,q2); */
  /* AimBr = MUL1(AimBr,q2); */
    
  /* ArpBi = MUL1(ArpBi,q2); */
  /* ArmBi = MUL1(ArmBi,q2); */
  
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
  
  
  X5 = SUB(ArmBi, AipBr);
  Y5 = ADD(AipBr, ArmBi);
  X7 = ADD(ArpBi, AimBr);
  Y7 = SUB(AimBr, ArpBi);
 
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

  //w = exp(2 pi i / 16)

  c = BROADCAST(1556715293UL);
  s = BROADCAST(978592373UL);

  AipBr = ADD(Ai, Br);				
  AimBr = SUB(Ai, Br);				
  ArpBi = ADD(Ar, Bi);				
  ArmBi = SUB(Ar, Bi);				
  t0 = mul1(ArmBi,c);				
  t1 = mul1(AipBr,s);				
  X9 = SUB(t0, t1);				
 			
  t0 = mul1(AipBr,c);				
  t1 = mul1(ArmBi,s);				
  Y9 = ADD(t0, t1);				
 		
  t0 = mul1(ArpBi, c);				
  t1 = mul1(AimBr, s);				
  X15 = ADD(t0, t1);				
  		
  t0 = mul1(AimBr,c);				
  t1 = mul1(ArpBi,s);				
  Y15 = SUB(t0, t1);				
		

  //w = exp(6 pi i / 16) = (s_16, c_16)

  AipBr = ADD(Ei, Fr);				
  AimBr = SUB(Ei, Fr);				
  ArpBi = ADD(Er, Fi);				
  ArmBi = SUB(Er, Fi);				
  t0 = mul1(ArmBi,s);				
  t1 = mul1(AipBr,c);				
  X13 = SUB(t0, t1);				
  
  t0 = mul1(AipBr,s);				
  t1 = mul1(ArmBi,c);				
  Y13 = ADD(t0, t1);				
  
  t0 = mul1(ArpBi,s);				
  t1 = mul1(AimBr,c);				
  X11 = ADD(t0, t1);				
  
  t0 = mul1(AimBr,s);				
  t1 = mul1(ArpBi,c);				
  Y11 = SUB(t0, t1);				
    

  //C4

  c4_reg_to_mem(X0, Y0, X1, Y1, X2, Y2, X3, Y3, a, a+1*V, a+2*V, a+3*V, a+4*V, a+5*V, a+6*V, a+7*V);

  //  C4_REG_MEM(X0, Y0, X1, Y1, X2, Y2, X3, Y3, a[0], a[V], a[2*V], a[3*V], a[4*V], a[5*V], a[6*V], a[7*V]);


  //C2 - 6 times

  c2_reg_to_mem(X4, Y4, X5, Y5, a+8*V, a+9*V, a+10*V, a+11*V);
 
  c2_reg_to_mem(X6, Y6, X7, Y7, a+12*V, a+13*V, a+14*V, a+15*V);

  
  c2_reg_to_mem(X8, Y8, X9, Y9, a+16*V, a+17*V, a+18*V, a+19*V);
  
  c2_reg_to_mem(X10, Y10, X11, Y11, a+20*V, a+21*V, a+22*V, a+23*V);
   
  c2_reg_to_mem(X12, Y12, X13, Y13, a+24*V, a+25*V, a+26*V, a+27*V);
    
  c2_reg_to_mem(X14, Y14, X15, Y15, a+28*V, a+29*V, a+30*V, a+31*V);
   
}
