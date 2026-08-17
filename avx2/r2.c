#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void r1024(register integer *a)
{
  rpass_srl3(a, d512, d1024, 128);

  r256(a);
  cs128(a + 256);
  cs128(a + 512);
  cs128(a + 768);
}
