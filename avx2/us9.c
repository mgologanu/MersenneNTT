#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void us131072(register integer *a)
{
  us65536(a);
  us32768(a + 131072);
  us32768(a + 196608);
  uspass_sr(a, d131072, 16384);
}
