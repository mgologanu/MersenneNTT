#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void us32768(register integer *a)
{
  us16384(a);
  us8192(a + 32768);
  us8192(a + 49152);
  uspass_sr(a, d32768, 4096);
}
