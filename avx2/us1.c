#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void us512(register integer *a)
{
  us128(a);
  us64(a + 256);
  us64(a + 384);
  us64(a + 512);
  us64(a + 640);
  us64(a + 768);
  us64(a + 896);
  uspass_srl3(a, d256, d512, 64);
}
