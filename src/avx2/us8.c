#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void us65536(register integer *a)
{
  us32768(a);
  us16384(a + 65536);
  us16384(a + 98304);
  uspass_sr(a, d65536, 8192);
}
