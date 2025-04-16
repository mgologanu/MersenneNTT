#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void v262144(register integer *a)
{
  v65536(a);
  us32768(a + 65536);
  us32768(a + 131072);
  us32768(a + 196608);

  vpass_srl3(a, d131072, d262144, 32768);
}
