#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void cs32768(register integer *a)
{
  cspass_sr(a, d32768, 4096);
  cs16384(a);
  cs8192(a + 32768);
  cs8192(a + 49152);
}
