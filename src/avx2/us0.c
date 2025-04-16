#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void  us16(register integer *a)
{
  u4_vert(a);

  uspass_rad4(a, rad8_d32, 2);
  
}

void  us32(register integer *a)
{
  u4_vert(a);   //Note that c4_vert also contains the transpose part

  u4_vert(a+32);

  uspass_rad8(a, rad8_d32, 4);

}


void  us64(register integer *a)
{
  us8_vert(a); //Note that c8_vert also contains the transpose part
  //  uts8_V4  (a);

  us8_vert(a+64);
  //  uts8_V4  (a+64);
  
  uspass_rad8(a, rad8_d64, 8);
}

void  us128(register integer *a)
{
  us16_vert(a);
  uts16_V4  (a);

  us16_vert(a+128);
  uts16_V4  (a+128);

  uspass_rad8(a, rad8_d128, 16);
}

void us256(register integer *a)
{
  us64(a);
  us32(a + 128);
  us32(a + 192);
  us32(a + 256);
  us32(a + 320);
  us32(a + 384);
  us32(a + 448);
  uspass_srl3(a, d128, d256, 32);
}
