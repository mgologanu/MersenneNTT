#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void r131072(register integer *a)
{
  rpass_srl3(a, d65536, d131072, 16384);

  r32768(a);
  cs16384(a + 32768);
  cs16384(a + 65536);
  cs16384(a + 98304);
}
