#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void v2048(register integer *a)
{
  v512(a);
  us256(a + 512);
  us256(a + 1024);
  us256(a + 1536);

  vpass_srl3(a, d1024, d2048, 256);
}
