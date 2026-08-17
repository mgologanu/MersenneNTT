#include "ntt.h"

#include "pre8_V4.h"
#include "i.h"


void cs1024(register integer *a)
{
  cspass_srl3(a, d512, d1024, 128);
  cs256(a);
  cs128(a + 512);
  cs128(a + 768);
  cs128(a + 1024);
  cs128(a + 1280);
  cs128(a + 1536);
  cs128(a + 1792);
}
