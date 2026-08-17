#include <stdalign.h>
#include <assert.h>
#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"
 
const _Alignas(32) integer d32768[] = {
#include "roots/p4_31_32768.txt"
};

static_assert(sizeof(d32768)/sizeof(integer) == 16384, "Wrong initializer list for d32768");
