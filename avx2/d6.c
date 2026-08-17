#include <stdalign.h>
#include <assert.h>
#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"
 
const _Alignas(32) integer d16384[] = {
#include "roots/p4_31_16384.txt"
};

static_assert(sizeof(d16384)/sizeof(integer) == 8192, "Wrong initializer list for d16384");
