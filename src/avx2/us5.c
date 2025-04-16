#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void us8192(register integer *a)
{
  us4096(a);
  us2048(a + 8192);
  us2048(a + 12288);
  uspass_sr(a, d8192, 1024);
}
