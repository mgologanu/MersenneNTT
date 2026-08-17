#include "ntt.h"

#include "pre8_V4.h"
#include "i.h"



void  cs16(register integer *a)
{

  /*
    the first half of rad8_32: w0, w, w2, w3

    (with n = 32/8 = 4)

    w0 = 1/2 of roots of order 2n = 8  and  
    w  = 1/4 of roots of order 4n = 16


    are precisely the roots needed for rad4_16:

    (with n = 16/8 = 2)
  
    w0 =  1/2 of roots of order 4n = 8
    w  =  1/4 of roots of order 8n = 16
 
  */

  
  cspass_rad4(a, rad8_d32, 2); 
  
  c4_vert(a);   //Note that c4_vert also contains the transpose part

}

void  cs32(register integer *a)
{
  cspass_rad8(a, rad8_d32, 4);

  c4_vert(a);   //Note that c4_vert also contains the transpose part

  c4_vert(a+32);

}


void  cs64(register integer *a)
{
  cspass_rad8(a, rad8_d64, 8);

  cs8_vert(a);
  
  cs8_vert(a+64);
}

void  cs128(register integer *a)
{
  cspass_rad8(a, rad8_d128, 16);

  ts16_V4  (a);
  cs16_vert(a);

  ts16_V4  (a+128);
  cs16_vert(a+128);
};



/* void  cs128(register integer *a) */
/* { */
/*   cspass_sr(a, d128, 16); */
/*   cs64(a); */
/*   cs32(a + 128); */
/*   cs32(a + 192); */
/* } */


void cs256(register integer *a)
{
  cspass_srl3(a, d128, d256, 32);
  cs64(a);
  cs32(a + 128);
  cs32(a + 192);
  cs32(a + 256);
  cs32(a + 320);
  cs32(a + 384);
  cs32(a + 448);
}

