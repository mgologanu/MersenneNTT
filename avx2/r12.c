#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void r1048576(register integer *a)
{
  rpass_srl3(a, d524288, d1048576, 131072);

  r262144(a);
  cs131072(a + 262144);
  cs131072(a + 524288);
  cs131072(a + 786432);
}
