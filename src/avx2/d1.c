#include <stdalign.h>
#include <assert.h>
#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"
 

const _Alignas(32) integer d512[] = {
#include "roots/p4_31_512.txt"
};
static_assert(sizeof(d512)/sizeof(integer) == 256, "Wrong initializer list for d256");

