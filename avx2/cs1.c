#include "ntt.h"

#include "pre8_V4.h"
#include "i.h"



void cs512(register integer *a)
{
  cspass_srl3(a, d256, d512, 64);
  cs128(a);
  cs64(a + 256);
  cs64(a + 384);
  cs64(a + 512);
  cs64(a + 640);
  cs64(a + 768);
  cs64(a + 896);
}
