#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void r8192(register integer *a)
{
  rpass_srl3(a, d4096, d8192, 1024);

  r2048(a);
  cs1024(a + 2048);
  cs1024(a + 4096);
  cs1024(a + 6144);
}
