#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void v8192(register integer *a)
{
  v2048(a);
  us1024(a + 2048);
  us1024(a + 4096);
  us1024(a + 6144);

  vpass_srl3(a, d4096, d8192, 1024);
}
