#include "ntt.h"
#include "pre8_V4.h"

#include "i.h"

/* 
   FFT8V - 8*V complex numbers

   
   V <= 8

   a: 8*V complex = 16*V values  interleaved V-packed format  - V real parts, V imag parts, etc.

*/
void us8_vert(integer *a)
{

  //n = N/8 = 1

  const integer_packed p_64 = BROADCAST(PRIME_64); 
  
  integer_packed x0, x1, x2, x3, x4, x5, x6, x7;
  integer_packed y0, y1, y2, y3, y4, y5, y6, y7;
  integer_packed Ar, Br, Ai, Bi, AipBr, AimBr, ArpBi, ArmBi;
  integer_packed Er, Ei, Fr, Fi;
  integer_packed ar,br,cr,dr,ai,bi,ci,di;
  integer_packed X0,Y0, Xm01, Ym01;
  integer_packed t0, t1, t2, t3, t4, t5;
  integer_packed s0, s1, s2, s3, s4, s5;
  /* 
        ---C8 -----  ---C8 -----  ---C8 ----- ---C8 ----- 
	x0 x1 y0 y1  x2 x3 y2 y3  x4 x5 y4 y5   x6 x7 y6 y7
    (V*) 0  1  2  3   4  5  6  7   8  9 10 11  12 13 14 15
     
     transpose(x0 x2 x4 x6) -> t0 t1 t2 t3 
     
     transpose(x1 x3 x5 x7) -> t4 t5 t6 t7
     
     save in place in increasing order of ti

     Same for y
      
  */

  //c = 1, s=0

  
  ArmBi = LOAD(a[8*V]); 
  AipBr = LOAD(a[9*V]);
  ArpBi = LOAD(a[14*V]);    
  AimBr = LOAD(a[15*V]);

  Ai = ADD(AipBr, AimBr);
  Br = SUB(AipBr, AimBr);

  Ar = ADD(ArpBi, ArmBi);
  Bi = SUB(ArpBi, ArmBi);
  
  //c = 1, s = 0
  
  ArmBi = LOAD(a[12*V]); 
  AipBr = LOAD(a[13*V]);
  ArpBi = LOAD(a[10*V]);    
  AimBr = LOAD(a[11*V]);
  

  Ei = ADD(AipBr, AimBr);
  Fr = SUB(AipBr, AimBr);

  Er = ADD(ArpBi, ArmBi);
  Fi = SUB(ArpBi, ArmBi);
  
 

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

  
  // c = 1, s = 0

  ArmBi = LOAD(a[4*V]); 
  AipBr = LOAD(a[5*V]);
  ArpBi = LOAD(a[6*V]);    
  AimBr = LOAD(a[7*V]);

  Ai = ADD(AipBr, AimBr);
  Br = SUB(AipBr, AimBr);

  Ar = ADD(ArpBi, ArmBi);
  Bi = SUB(ArpBi, ArmBi);

  // c = 1, s = 0

  Xm01 = LOAD(a[2*V]); //X 0+4+2+6-(1+5+3+7)
  Ym01 = LOAD(a[3*V]); //Y 0+4+2+6-(1+5+3+7)

  X0 = LOAD(a[0]);     //X 0+4+2+6+(1+5+3+7)
  Y0 = LOAD(a[V]);     //Y 0+4+2+6+(1+5+3+7)


  
  t0 = ADD(X0, Xm01); //X 0+4+2+6
  s0 = ADD(Y0, Ym01); //Y 0+4+2+6
  
  
  t1 = SUB(X0, Xm01); //X 1+5+3+7
  s1 = SUB(Y0, Ym01); //Y 1+5+3+7
  
  t2 = ADD(t0, Ar); //X 0+4
  t3 = SUB(t0, Ar); //X 2+6
  
  s2 = ADD(s0, Ai); //Y 0+4
  s3 = SUB(s0, Ai); //Y 2+6

  x0 = ADD(t2, ar);
  x4 = SUB(t2, ar);
  
  y0 = ADD(s2, ai);
  y4 = SUB(s2, ai);

  
  x2 = ADD(t3, br);
  x6 = SUB(t3, br);
  
  y2 = ADD(s3, bi);
  y6 = SUB(s3, bi);

  
  t2 = ADD(t1, Br); //X 1+5
  t3 = SUB(t1, Br); //X 3+7
  
  s2 = ADD(s1, Bi); //Y 1+5
  s3 = SUB(s1, Bi); //Y 3+7
  
  x1 = ADD(t2, cr);
  x5 = SUB(t2, cr);
  
  y1 = ADD(s2, ci);
  y5 = SUB(s2, ci);

  
   
  x3 = ADD(t3, dr);
  x7 = SUB(t3, dr);
  
  y3 = ADD(s3, di);
  y7 = SUB(s3, di);


  //TRANSPOSE4(x0,x1,x2,x3);

  transpose4(&x0,&x1,&x2,&x3);
      
  STORE0(a[0],    x0);
  STORE0(a[4*V],  x1);
  STORE0(a[8*V],  x2);
  STORE0(a[12*V], x3);

  //  TRANSPOSE4(x4,x5,x6,x7);
  transpose4(&x4,&x5,&x6,&x7);

  STORE0(a[1*V],  x4);
  STORE0(a[5*V],  x5);
  STORE0(a[9*V], x6);
  STORE0(a[13*V], x7);

  //  TRANSPOSE4(y0,y1,y2,y3);
  transpose4(&y0,&y1,&y2,&y3);
  
  STORE0(a[2*V],  y0);
  STORE0(a[6*V],  y1);
  STORE0(a[10*V], y2);
  STORE0(a[14*V], y3);

  //  TRANSPOSE4(y4,y5,y6,y7);
  transpose4(&y4,&y5,&y6,&y7);

  STORE0(a[3*V],  y4);
  STORE0(a[7*V],  y5);
  STORE0(a[11*V], y6);
  STORE0(a[15*V], y7);


}
