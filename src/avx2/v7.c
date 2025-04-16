#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void v32768(register integer *a)
{
  v8192(a);
  us4096(a + 8192);
  us4096(a + 16384);
  us4096(a + 24576);

  vpass_srl3(a, d16384, d32768, 4096);
}
