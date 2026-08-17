#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void v16384(register integer *a)
{
  v4096(a);
  us2048(a + 4096);
  us2048(a + 8192);
  us2048(a + 12288);

  vpass_srl3(a, d8192, d16384, 2048);
}
