#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

void cs1048576(register integer *a)
{
  cspass_sr(a, d1048576, 131072);
  cs524288(a);
  cs262144(a + 1048576);
  cs262144(a + 1572864);
}
