#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void r512(register integer *a)
{
  rpass_srl3(a, d256, d512, 64);

  r128(a);
  cs64(a + 128);
  cs64(a + 256);
  cs64(a + 384);
}
