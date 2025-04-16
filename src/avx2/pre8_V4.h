#include<immintrin.h>


#define V 4

#define V2 8

#define V4 16

#define V8 32

#define VPOW 2


#define PRIME   2147483647UL

//#define PRIME_64 137438953408UL

#define PRIME_64 274877906816UL

#define P_POW 31

#define SQ2_POW 15

#define ADD(a,b)     _mm256_add_epi64(a,b)
#define SUB(a,b)     _mm256_sub_epi64(a,b)
#define LOAD(a)      _mm256_load_si256((__m256i *) &a)
#define LOADU(a)     _mm256_loadu_si256((__m256i *) &a)
#define BROADCAST(a) _mm256_set1_epi64x(a)


#define STORE0(a,b)   _mm256_store_si256((__m256i *) &a,b)

#define STORE(a,b)  {				\
    red1(&b);				\
    _mm256_store_si256((__m256i *) &a,b);	\
  }


#define STORE1(a,b)   _mm256_store_si256((__m256i *) a,b)


#define MUL0(a, b)     ({				\
      _mul0_t1  = _mm256_mul_epu32(a, b);			\
      _mul0_t2  = _mm256_and_si256(_mul0_t1, p);			\
      _mul0_t3  = _mm256_srli_epi64(_mul0_t1, P_POW);		\
      retval = _mm256_add_epi64(_mul0_t2, _mul0_t3);		\
      retval;						\
    })

/*
  a = in interval [-k*p, k*p] with k = 5 + 2^stages < 64

  Add 64*p to a so that it is positive

  Do a reduction to [0 p+64] aproximative, inside [0 2p]

  Do the multiplication with b assumed in [0 2p]

  Do a single reduction to get back to [0 5p]
  
 */


#define MUL1(a, b)     ({				\
      _mul1_t3  = _mm256_add_epi64(a, p_64);			\
      _mul1_t1  = _mm256_and_si256(_mul1_t3, p);			\
      _mul1_t2  = _mm256_srli_epi64(_mul1_t3, P_POW);		\
      _mul1_t3  = _mm256_add_epi64(_mul1_t1, _mul1_t2);		\
      _mul1_t1  = _mm256_mul_epu32(_mul1_t3, b);			\
      _mul1_t2  = _mm256_and_si256(_mul1_t1, p);			\
      _mul1_t3  = _mm256_srli_epi64(_mul1_t1, P_POW);		\
      retval = _mm256_add_epi64(_mul1_t2, _mul1_t3);		\
      retval;						\
    })



/* #define MUL2(a, b)     ({					\ */
/*       _t1  = _mm256_mul_epu32(a, b);				\ */
/*       _t2  = _mm256_and_si256(_t1, p);				\ */
/*       _t3  = _mm256_srli_epi64(_t1, P_POW);			\ */
/*       _t1  = _mm256_add_epi64(_t2, _t3);			\ */
/*       _t2  = _mm256_and_si256(_t1, p);				\ */
/*       _t3  = _mm256_srli_epi64(_t1, P_POW);			\ */
/*       retval = _mm256_add_epi64(_t2, _t3);			\ */
/*       retval;							\ */
/*     }) */


/* #define MUL3(a, b)     ({					\ */
/*       _t1  = _mm256_mul_epu32(a, b);				\ */
/*       _t2  = _mm256_and_si256(_t1, p);				\ */
/*       _t3  = _mm256_srli_epi64(_t1, P_POW);			\ */
/*       _t1  = _mm256_add_epi64(_t2, _t3);			\ */
/*       _t2  = _mm256_and_si256(_t1, p);				\ */
/*       _t3  = _mm256_srli_epi64(_t1, P_POW);			\ */
/*       _t3  = _mm256_add_epi64(_t2, _t3);			\ */
/*       _t1 =  _mm256_sub_epi64(_t3, p);				\ */
/*       _tk = _mm256_cmpgt_epi64 (p, _t1);			\ */
/*       retval =  _mm256_blendv_epi8(_t3, _t1, _tk);		\ */
/*       retval;							\ */
/*     }) */

/*
  T = in interval [-k*p, k*p] with k = 5 + 2^stages < 64

  Add 64*p to a so that it is positive

  Do a reduction to [0 p+64] aproximative, inside [0 2p]
*/

#define RED1(T)     ({							\
      _red1_t3  = _mm256_add_epi64(       T,    p_64);			\
      _red1_t1  = _mm256_and_si256( _red1_t3,   p);			\
      _red1_t2  = _mm256_srli_epi64(_red1_t3,   P_POW);			\
      retval    = _mm256_add_epi64( _red1_t1,  _red1_t2);		\
      retval;								\
    })


#define RED2(T)     ({							\
      _red2_t3  = _mm256_add_epi64(T, p_64);				\
      _red2_t1  = _mm256_and_si256(_red2_t3, p);			\
      _red2_t2  = _mm256_srli_epi64(_red2_t3, P_POW);			\
      _red2_t3  = _mm256_add_epi64(_red2_t1, _red2_t2);			\
      _red2_t1  = _mm256_and_si256(_red2_t3, p);			\
      _red2_t2  = _mm256_srli_epi64(_red2_t3, P_POW);			\
      retval = _mm256_add_epi64(_red2_t1, _red2_t2);			\
      retval;								\
    })


#define RED3(T)     ({							\
      _red3_t3  = _mm256_add_epi64(T, p_64);				\
      _red3_t1  = _mm256_and_si256(_red3_t3, p);			\
      _red3_t2  = _mm256_srli_epi64(_red3_t3, P_POW);			\
      _red3_t3  = _mm256_add_epi64(_red3_t1, _red3_t2);			\
      _red3_t1  = _mm256_and_si256(_red3_t3, p);			\
      _red3_t2  = _mm256_srli_epi64(_red3_t3, P_POW);			\
      _red3_t3  =  _mm256_add_epi64(_red3_t1, _red3_t2);		\
      _red3_t1  =  _mm256_sub_epi64(_red3_t3, p);			\
      _red3_tk  = _mm256_cmpgt_epi64 (p, _red3_t3);			\
      retval =  _mm256_blendv_epi8(_red3_t1, _red3_t3, _red3_tk);	\
      retval;								\
    })




#define TWISTED(Ar,Ai,Br,Bi,c,s,S0,S1,S2,S3) {		\
    _AipBr = ADD(Ai, Br);				\
    _AimBr = SUB(Ai, Br);				\
    _ArpBi = ADD(Ar, Bi);				\
    _ArmBi = SUB(Ar, Bi);				\
    _t0 = MUL1(_ArmBi, c);				\
    _t1 = MUL1(_AipBr, s);				\
    _t2 = SUB(_t0, _t1);				\
    STORE(S0, _t2);					\
    _t0 = MUL1(_AipBr, c);				\
    _t1 = MUL1(_ArmBi, s);				\
    _t2 = ADD(_t0, _t1);				\
    STORE(S1, _t2);					\
    _t0 = MUL1(_ArpBi, c);				\
    _t1 = MUL1(_AimBr, s);				\
    _t2 = ADD(_t0, _t1);				\
    STORE(S2, _t2);					\
    _t0 = MUL1(_AimBr, c);				\
    _t1 = MUL1(_ArpBi, s);				\
    _t2 = SUB(_t0, _t1);				\
    STORE(S3, _t2);					\
  }




/*                    x0   x1   x2   x3    */
/* offset + stride * ( 0    1    2    3  ) */

#define TRANSPOSE4_FROM_MEM(x0,x1,x2,x3,a,offset,stride ) {		\
    _t0 = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_load_si128((__m128i *) &a[offset+0*stride+0])), _mm_load_si128((__m128i *) &a[offset+2*stride+0]), 1); \
    _t1 = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_load_si128((__m128i *) &a[offset+1*stride+0])), _mm_load_si128((__m128i *) &a[offset+3*stride+0]), 1); \
    _t2 = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_load_si128((__m128i *) &a[offset+0*stride+2])), _mm_load_si128((__m128i *) &a[offset+2*stride+2]), 1); \
    _t3 = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_load_si128((__m128i *) &a[offset+1*stride+2])), _mm_load_si128((__m128i *) &a[offset+3*stride+2]), 1); \
    x0 = _mm256_unpacklo_epi64(_t0,_t1);				\
    x1 = _mm256_unpackhi_epi64(_t0,_t1);				\
    x2 = _mm256_unpacklo_epi64(_t2,_t3);					\
    x3 = _mm256_unpackhi_epi64(_t2,_t3);					\
 }

#define PERM(a,b,p)  _mm256_permute2f128_si256(a,b,p)
#define SHUF(a,b,i)  (__m256i ) _mm256_shuffle_pd((__m256d) a,(__m256d) b,i)

#define P1_F 0x20
#define P2_F 0x31
#define I1_F 0x0
#define I2_F 0xF


#define TRANSPOSE4(x0,x1,x2,x3) {				   \
    _s0 = SHUF(x0,x1,I1_F) ;                                       \
    _s1 = SHUF(x0,x1,I2_F) ;                                       \
    _s2 = SHUF(x2,x3,I1_F) ;                                       \
    _s3 = SHUF(x2,x3,I2_F) ;					    \
    x0 = PERM(_s0,_s2,P1_F) ;                                       \
    x2 = PERM(_s0,_s2,P2_F) ;                                       \
    x1 = PERM(_s1,_s3,P1_F) ;                                       \
    x3 = PERM(_s1,_s3,P2_F) ;                                       \
  }


#define C4_REG_MEM(x0, y0, x1, y1, x2, y2, x3, y3, X0, Y0, X1, Y1, X2, Y2, X3, Y3) { \
    _t02 = ADD(x0, x2);          /* x 0+2 */				\
    _Ar  = SUB(x0, x2);          /* x 0-2 */				\
    _t13 = ADD(x1, x3);          /* x 1+3 */				\
    _Br  = SUB(x1, x3);          /* x 1-3 */				\
    _s02 = ADD(y0, y2);          /* y 0+2 */				\
    _Ai  = SUB(y0, y2);          /* y 0-2 */				\
    _s13 = ADD(y1, y3);          /* y 1+3 */				\
    _Bi  = SUB(y1, y3);          /* y 1-3 */				\
    _t0  = ADD(_t02, _t13);         /* x 0+2+1+3 */			\
    _t1  = ADD(_s02, _s13);         /* y 0+2+1+3 */			\
    _t2  = SUB(_t02, _t13);         /* x 0+2-1-3 */			\
    _t3  = SUB(_s02, _s13);         /* y 0+2-1-3 */			\
    _AipBr = ADD(_Ai, _Br);        /* y 0-2 + x 1-3 */			\
    _AimBr = SUB(_Ai, _Br);        /* y 0-2 - x 1-3 */			\
    _ArpBi = ADD(_Ar, _Bi);        /* x 0-2 + y 1-3 */			\
    _ArmBi = SUB(_Ar, _Bi);        /* x 0-2 - y 1-3 */			\
    STORE(X0, _t0);							\
    STORE(Y0, _t1);							\
    STORE(X1, _t2);							\
    STORE(Y1, _t3);							\
    STORE(X2, _ArmBi);							\
    STORE(Y2, _AipBr);							\
    STORE(X3, _ArpBi);							\
    STORE(Y3, _AimBr);							\
  }



#define C2_REG_MEM(x0, y0, x1, y1, X0, Y0, X1, Y1) { \
    _t0 = ADD(x0, x1);				     \
    _t2 = SUB(x0, x1);				     \
    _t1 = ADD(y0, y1);				     \
    _t3 = SUB(y0, y1);				     \
    STORE(X0, _t0);				     \
    STORE(Y0, _t1);				     \
    STORE(X1, _t2);				     \
    STORE(Y1, _t3);				     \
  } 
