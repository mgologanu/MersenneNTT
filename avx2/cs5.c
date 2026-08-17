#include "ntt.h"

#include "pre8_V4.h"
#include "i.h"



/* void cs8192(register integer *a) */
/* { */
/*   cspass_srl3(a, d4096, d8192, 1024); */
/*   cs2048(a); */
/*   cs1024(a + 4096); */
/*   cs1024(a + 6144); */
/*   cs1024(a + 8192); */
/*   cs1024(a + 10240); */
/*   cs1024(a + 12288); */
/*   cs1024(a + 14336); */
/* } */

void cs8192(register integer *a)
{
  cspass_sr(a, d8192, 1024);
  cs4096(a);
  cs2048(a + 8192);
  cs2048(a + 12288);
}
