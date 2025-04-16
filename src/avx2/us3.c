#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void us2048(register integer *a)
{
  us512(a);
  us256(a + 1024);
  us256(a + 1536);
  us256(a + 2048);
  us256(a + 2560);
  us256(a + 3072);
  us256(a + 3584);
  uspass_srl3(a, d1024, d2048, 256);
}
