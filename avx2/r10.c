#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void r262144(register integer *a)
{
  rpass_srl3(a, d131072, d262144, 32768);

  r65536(a);
  cs32768(a + 65536);
  cs32768(a + 131072);
  cs32768(a + 196608);
}
