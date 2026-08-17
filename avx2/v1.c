#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"


void v512(register integer *a)
{
  v128(a);
  us64(a + 128);
  us64(a + 256);
  us64(a + 384);

  vpass_srl3(a, d256, d512, 64);
}
