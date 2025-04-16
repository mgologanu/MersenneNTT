#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void us262144(register integer *a)
{
  us131072(a);
  us65536(a + 262144);
  us65536(a + 393216);
  uspass_sr(a, d262144, 32768);
}
