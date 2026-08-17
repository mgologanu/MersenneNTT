#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void r32768(register integer *a)
{
  rpass_srl3(a, d16384, d32768, 4096);

  r8192(a);
  cs4096(a + 8192);
  cs4096(a + 16384);
  cs4096(a + 24576);
}
