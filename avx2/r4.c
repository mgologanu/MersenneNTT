#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void r4096(register integer *a)
{
  rpass_srl3(a, d2048, d4096, 512);

  r1024(a);
  cs512(a + 1024);
  cs512(a + 2048);
  cs512(a + 3072);
}
