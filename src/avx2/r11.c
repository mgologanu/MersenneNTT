#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void r524288(register integer *a)
{
  rpass_srl3(a, d262144, d524288, 65536);

  r131072(a);
  cs65536(a + 131072);
  cs65536(a + 262144);
  cs65536(a + 393216);
}
