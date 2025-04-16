#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void us1024(register integer *a)
{
  us256(a);
  us128(a + 512);
  us128(a + 768);
  us128(a + 1024);
  us128(a + 1280);
  us128(a + 1536);
  us128(a + 1792);
  uspass_srl3(a, d512, d1024, 128);
}
