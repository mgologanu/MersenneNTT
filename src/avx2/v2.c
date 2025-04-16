#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void v1024(register integer *a)
{
  v256(a);
  us128(a + 256);
  us128(a + 512);
  us128(a + 768);

  vpass_srl3(a, d512, d1024, 128);
}
