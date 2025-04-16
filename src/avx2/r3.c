#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void r2048(register integer *a)
{
  rpass_srl3(a, d1024, d2048, 256);

  r512(a);
  cs256(a + 512);
  cs256(a + 1024);
  cs256(a + 1536);
}
