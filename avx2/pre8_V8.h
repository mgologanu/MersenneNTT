#include<immintrin.h>

typedef int64_t integer;
typedef __m512i integer_packed;

#define V 8

#define V2 16

#define VPOW 3

#define PRIME   2147483647L

#define PM1_HALF  1073741823L


#define ADD(a,b)     _mm512_add_epi64(a,b)
#define SUB(a,b)     _mm512_sub_epi64(a,b)
#define LOAD(a)      _mm512_load_epi64(&a)
#define LOADU(a)     _mm512_loadu_epi64(&a)
#define STORE(a,b)   _mm512_store_epi64(&a,b)
#define BROADCAST(a) _mm512_set1_epi64(a)


#define MUL1(a, b)     ({				\
      t1  = _mm512_mullo_epi64(a, b);			\
      t2  = _mm512_and_epi64(t1, p);			\
      t3  = _mm512_srai_epi64(t1, 31);			\
      retval = _mm512_add_epi64(t2, t3);		\
      retval;						\
    })

#define MUL2(ab, a, b)     {				\
    t1  = _mm512_mullo_epi64(a, b);			\
    t2  = _mm512_and_epi64(t1, p);			\
    t3  = _mm512_srai_epi64(t1, 31);			\
    t1  = _mm512_add_epi64(t2, t3);			\
    t2  = _mm512_and_epi64(t1, p);			\
    t3  = _mm512_srai_epi64(t1, 31);			\
    ab = _mm512_add_epi64(t2, t3);			\
  }

#define MUL3(ab, a, b)     {					\
    t1   = _mm512_mullo_epi64(a, b);				\
    t2  = _mm512_and_epi64(t1, p);				\
    t3  = _mm512_srai_epi64(t1, 31);				\
    t1  = _mm512_add_epi64(t2, t3);				\
    t2  = _mm512_and_epi64(t1, p);				\
    t3  = _mm512_srai_epi64(t1, 31);				\
    t1  = _mm512_add_epi64(t2, t3);				\
    t2  = _mm512_sub_epi64(t1, p);				\
    tk  = _mm512_cmp_epi64_mask (pm1_half, t1, _MM_CMPINT_LT);	\
    ab  =  _mm512_mask_blend_epi64 (tk, t1, t2);		\
  }

#define RED1(T)     {				\
    t1 = _mm512_and_epi64(T, p);		\
    t2 = _mm512_srai_epi64(T, 31);		\
    T  = _mm512_add_epi64(t1, t2);		\
  }

#define RED2(T)     {				\
    t1 = _mm512_and_epi64(T, p);		\
    t2 = _mm512_srai_epi64(T, 31);		\
    T  = _mm512_add_epi64(t1, t2);		\
    t1 = _mm512_and_epi64(T, p);		\
    t2 = _mm512_srai_epi64(T, 31);		\
    T  = _mm512_add_epi64(t1, t2);		\
  }

#define RED2_FINAL(T)     {						\
    t1 = _mm512_and_epi64(T, p);				\
    t2 = _mm512_srai_epi64(T, 31);				\
    T  = _mm512_add_epi64(t1, t2);				\
    t2 = _mm512_and_epi64(T, p);				\
    t3 = _mm512_srai_epi64(T, 31);				\
    t1  = _mm512_add_epi64(t2, t3);				\
    t2  = _mm512_sub_epi64(t1, p);				\
    tk  = _mm512_cmp_epi64_mask (pm1_half, t1, _MM_CMPINT_LT);	\
    T  =  _mm512_mask_blend_epi64 (tk, t1, t2);			\
  }



