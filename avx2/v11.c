#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void v524288(register integer *a)
{
  v131072(a);
  us65536(a + 131072);
  us65536(a + 262144);
  us65536(a + 393216);

  vpass_srl3(a, d262144, d524288, 65536);
}
