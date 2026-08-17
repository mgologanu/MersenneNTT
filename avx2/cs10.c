#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void cs262144(register integer *a)
{
  cspass_sr(a, d262144, 32768);
  cs131072(a);
  cs65536(a + 262144);
  cs65536(a + 393216);
}
