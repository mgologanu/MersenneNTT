#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void r65536(register integer *a)
{
  rpass_srl3(a, d32768, d65536, 8192);

  r16384(a);
  cs8192(a + 16384);
  cs8192(a + 32768);
  cs8192(a + 49152);
}
