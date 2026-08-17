#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void cs524288(register integer *a)
{
  cspass_sr(a, d524288, 65536);
  cs262144(a);
  cs131072(a + 524288);
  cs131072(a + 786432);
}
