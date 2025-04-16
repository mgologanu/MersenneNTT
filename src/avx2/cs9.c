#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void cs131072(register integer *a)
{
  cspass_sr(a, d131072, 16384);
  cs65536(a);
  cs32768(a + 131072);
  cs32768(a + 196608);
}
