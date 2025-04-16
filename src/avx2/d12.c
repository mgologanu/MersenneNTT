#include <stdalign.h>
#include <assert.h>
#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"
 
const _Alignas(32) integer d1048576[] = {
#include "roots/p4_31_1048576.txt"
};

static_assert(sizeof(d1048576)/sizeof(integer) == 524288, "Wrong initializer list for d1048576");
