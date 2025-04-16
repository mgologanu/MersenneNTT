#include "ntt.h"
#include "pre8_V4.h"

#include "i.h"

void u4_vert(integer *a)
{

  integer_packed x0, y0, x1, y1, x2, y2, x3, y3;


  //  U4_REG_MEM(a[0*V], a[1*V], a[2*V], a[3*V], a[4*V], a[5*V], a[6*V], a[7*V], x0, y0, x1, y1, x2, y2, x3, y3);

  u4_reg_to_mem(a+0*V, a+1*V, a+2*V, a+3*V, a+4*V, a+5*V, a+6*V, a+7*V, &x0, &y0, &x1, &y1, &x2, &y2, &x3, &y3);
  

  //  TRANSPOSE4(x0,x1,x2,x3);
  //  TRANSPOSE4(y0,y1,y2,y3);

  transpose4(&x0,&x1,&x2,&x3);
  transpose4(&y0,&y1,&y2,&y3);

  

  STORE0(a[0*V], x0);		
  STORE0(a[1*V], y0);		
  STORE0(a[2*V], x1);		
  STORE0(a[3*V], y1);		
  STORE0(a[4*V], x2);		
  STORE0(a[5*V], y2);		
  STORE0(a[6*V], x3);		
  STORE0(a[7*V], y3);		
  
}
