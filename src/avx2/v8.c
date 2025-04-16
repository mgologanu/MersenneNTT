#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void v65536(register integer *a)
{
  v16384(a);
  us8192(a + 16384);
  us8192(a + 32768);
  us8192(a + 49152);

  vpass_srl3(a, d32768, d65536, 8192);
}
