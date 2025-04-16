#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void v1048576(register integer *a)
{
  v262144(a);
  us131072(a + 262144);
  us131072(a + 524288);
  us131072(a + 786432);

  vpass_srl3(a, d524288, d1048576, 131072);
}
