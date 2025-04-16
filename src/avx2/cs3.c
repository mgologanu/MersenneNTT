#include "ntt.h"

#include "pre8_V4.h"
#include "i.h"



void cs2048(register integer *a)
{
  cspass_srl3(a, d1024, d2048, 256);
  cs512(a);
  cs256(a + 1024);
  cs256(a + 1536);
  cs256(a + 2048);
  cs256(a + 2560);
  cs256(a + 3072);
  cs256(a + 3584);
}

