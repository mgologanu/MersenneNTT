#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void us4096(register integer *a)
{
  us2048(a);
  us1024(a + 4096);
  us1024(a + 6144);
  uspass_sr(a, d4096, 512);
}
