#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void r16384(register integer *a)
{
  rpass_srl3(a, d8192, d16384, 2048);

  r4096(a);
  cs2048(a + 4096);
  cs2048(a + 8192);
  cs2048(a + 12288);
}
