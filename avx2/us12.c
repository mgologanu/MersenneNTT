#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void us1048576(register integer *a)
{
  us524288(a);
  us262144(a + 1048576);
  us262144(a + 1572864);
  uspass_sr(a, d1048576, 131072);
}
