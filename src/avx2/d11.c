#include <stdalign.h>
#include <assert.h>
#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"
 
const _Alignas(32) integer d524288[] = {
#include "roots/p4_31_524288.txt"
};

static_assert(sizeof(d524288)/sizeof(integer) == 262144, "Wrong initializer list for d524288");
