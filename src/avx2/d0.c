#include <stdalign.h>
#include <assert.h>

#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"



const _Alignas(64) integer rad8_d32[] = {
#include "roots/p4_31_rad8_32.txt"
};

static_assert(sizeof(rad8_d32)/sizeof(integer) == 32, "Wrong initializer list for rad8_32");

const _Alignas(64) integer rad8_d64[] = {
#include "roots/p4_31_rad8_64.txt"
};

static_assert(sizeof(rad8_d64)/sizeof(integer) == 64, "Wrong initializer list for rad8_64");

const _Alignas(64) integer rad8_d128[] = {
#include "roots/p4_31_rad8_128.txt"
};
static_assert(sizeof(rad8_d128)/sizeof(integer) == 128, "Wrong initializer list for rad8_128");



/* const _Alignas(64) integer d16[] = { */
/* #include "roots/p4_31_16.txt" */
/* }; */

/* static_assert(sizeof(d128)/sizeof(integer) == 8, "Wrong initializer list for d16"); */


const _Alignas(64) integer d16[] = {
#include "roots/p4_31_16.txt"
};
static_assert(sizeof(d16)/sizeof(integer) == 8, "Wrong initializer list for d16");

const _Alignas(64) integer d32[] = {
#include "roots/p4_31_32.txt"
};
static_assert(sizeof(d32)/sizeof(integer) == 16, "Wrong initializer list for d32");

const _Alignas(64) integer d64[] = {
#include "roots/p4_31_64.txt"
};
static_assert(sizeof(d64)/sizeof(integer) == 32, "Wrong initializer list for d64");

const _Alignas(64) integer d128[] = {
#include "roots/p4_31_128.txt"
};
static_assert(sizeof(d128)/sizeof(integer) == 64, "Wrong initializer list for d128");

const _Alignas(64) integer d256[] = {
#include "roots/p4_31_256.txt"
};

static_assert(sizeof(d256)/sizeof(integer) == 128, "Wrong initializer list for d256");

