#include <stdalign.h>
#include <assert.h>
#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"
 
const _Alignas(32) integer d262144[] = {
#include "roots/p4_31_262144.txt"
};

static_assert(sizeof(d262144)/sizeof(integer) == 131072, "Wrong initializer list for d262144");
