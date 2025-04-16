#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void us16384(register integer *a)
{
  us8192(a);
  us4096(a + 16384);
  us4096(a + 24576);
  uspass_sr(a, d16384, 2048);
}
