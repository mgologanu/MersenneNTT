#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void cs65536(register integer *a)
{
  cspass_sr(a, d65536, 8192);
  cs32768(a);
  cs16384(a + 65536);
  cs16384(a + 98304);
}
