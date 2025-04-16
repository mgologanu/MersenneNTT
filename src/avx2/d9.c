#include <stdalign.h>
#include <assert.h>
#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"
 
const _Alignas(32) integer d131072[] = {
#include "roots/p4_31_131072.txt"
};

static_assert(sizeof(d131072)/sizeof(integer) == 65536, "Wrong initializer list for d131072");
