#ifndef _NTT_H
#define _NTT_H       


#include<immintrin.h>
#include <stdint.h>


typedef uint64_t integer;
typedef __m256i integer_packed;



extern void cs16(integer *a);
extern void cs32(integer *a);
extern void cs64(integer *a);
extern void cs128(integer *a);
extern void cs256(integer *a);
extern void cs512(integer *a);
extern void cs1024(integer *a);
extern void cs2048(integer *a);
extern void cs4096(integer *a);
extern void cs8192(integer *a);
extern void cs16384(integer *a);
extern void cs32768(integer *a);
extern void cs65536(integer *a);
extern void cs131072(integer *a);
extern void cs262144(integer *a);
extern void cs524288(integer *a);
extern void cs1048576(integer *a);



extern void us16(integer *a);
extern void us32(integer *a);
extern void us64(integer *a);
extern void us128(integer *a);
extern void us256(integer *a);
extern void us512(integer *a);
extern void us1024(integer *a);
extern void us2048(integer *a);
extern void us4096(integer *a);
extern void us8192(integer *a);
extern void us16384(integer *a);
extern void us32768(integer *a);
extern void us65536(integer *a);
extern void us131072(integer *a);
extern void us262144(integer *a);
extern void us524288(integer *a);
extern void us1048576(integer *a);


extern void r16(integer *a);
extern void r32(integer *a);
extern void r64(integer *a);
extern void r128(integer *a);
extern void r256(integer *a);
extern void r512(integer *a);
extern void r1024(integer *a);
extern void r2048(integer *a);
extern void r4096(integer *a);
extern void r8192(integer *a);
extern void r16384(integer *a);
extern void r32768(integer *a);
extern void r65536(integer *a);
extern void r131072(integer *a);
extern void r262144(integer *a);
extern void r524288(integer *a);
extern void r1048576(integer *a);


extern void v16(integer *a);
extern void v32(integer *a);
extern void v64(integer *a);
extern void v128(integer *a);
extern void v256(integer *a);
extern void v512(integer *a);
extern void v1024(integer *a);
extern void v2048(integer *a);
extern void v4096(integer *a);
extern void v8192(integer *a);
extern void v16384(integer *a);
extern void v32768(integer *a);
extern void v65536(integer *a);
extern void v131072(integer *a);
extern void v262144(integer *a);
extern void v524288(integer *a);
extern void v1048576(integer *a);


extern void mulc(integer *a, integer *b, size_t n);
extern void mulr(integer *a, integer *b, size_t n);

extern void scalec(integer *a, size_t n);
extern void scaler(integer *a, size_t n);
#endif /* ntt.h */                
