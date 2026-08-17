#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void us524288(register integer *a)
{
  us262144(a);
  us131072(a + 524288);
  us131072(a + 786432);
  uspass_sr(a, d524288, 65536);
}
