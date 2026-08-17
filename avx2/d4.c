#include <stdalign.h>
#include <assert.h>
#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"
 
const _Alignas(32) integer d4096[] = {
#include "roots/p4_31_4096.txt"
};

static_assert(sizeof(d4096)/sizeof(integer) == 2048, "Wrong initializer list for d4096");
