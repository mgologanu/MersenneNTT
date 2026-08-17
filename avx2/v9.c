#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void v131072(register integer *a)
{
  v32768(a);
  us16384(a + 32768);
  us16384(a + 65536);
  us16384(a + 98304);

  vpass_srl3(a, d65536, d131072, 16384);
}
