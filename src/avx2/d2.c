#include <stdalign.h>
#include <assert.h>
#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"
 
const _Alignas(32) integer d1024[] = {
#include "roots/p4_31_1024.txt"
};

static_assert(sizeof(d1024)/sizeof(integer) == 512, "Wrong initializer list for d1024");
