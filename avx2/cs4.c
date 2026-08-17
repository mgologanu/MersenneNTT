#include "ntt.h"

#include "pre8_V4.h"
#include "i.h"



/* void cs4096(register integer *a) */
/* { */
/*   cspass_srl3(a, d2048, d4096, 512); */
/*   cs1024(a); */
/*   cs512(a + 2048); */
/*   cs512(a + 3072); */
/*   cs512(a + 4096); */
/*   cs512(a + 5120); */
/*   cs512(a + 6144); */
/*   cs512(a + 7168); */
/* } */


void cs4096(register integer *a)
{
  cspass_sr(a, d4096, 512);
  cs2048(a);
  cs1024(a + 4096);
  cs1024(a + 6144);
}

