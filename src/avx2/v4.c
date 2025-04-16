#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void v4096(register integer *a)
{
  v1024(a);
  us512(a + 1024);
  us512(a + 2048);
  us512(a + 3072);

  vpass_srl3(a, d2048, d4096, 512);
}
