#ifndef _MRSN_AVX2_INTERNAL_H
#define _MRSN_AVX2_INTERNAL_H


#define MRSN_PRIME 2147483647

#define MRSN_PRIME_POW 31

static inline __attribute__((always_inline))
void add_sub(__m256i * a, __m256i * b)
{

  __m256i t1, t2, minus_b;

  const __m256i p = _mm256_set1_epi32(MRSN_PRIME);
  


  minus_b = _mm256_xor_si256(*b, p); //~b

  t1 = _mm256_add_epi32(*a, *b);     //a+b
  t2 = _mm256_add_epi32(*a, minus_b);//a-b

  *a = _mm256_and_si256(t1, p);//a+b - 30:0 bits
  *b = _mm256_and_si256(t2, p);//a-b - 30:0 bits

  t1 = _mm256_srli_epi32(t1, MRSN_PRIME_POW);//a+b >> 31
  t2 = _mm256_srli_epi32(t2, MRSN_PRIME_POW);//a+b >> 31

  
  *a = _mm256_add_epi32(t1, *a);
  *b = _mm256_add_epi32(t2, *b);

  //Overflow? Only if a+b = 0xFFFFFFFF = 2*p+1.
  //But a<=p, b<=p, overflow here never happens
    
}

static inline __attribute__((always_inline))
__m256i add(__m256i a, __m256i b)
{

  __m256i t1, t2, t3;

  const __m256i p = _mm256_set1_epi32(MRSN_PRIME);
  
  t3 = _mm256_add_epi32(a, b);
  
  t1 = _mm256_and_si256(t3, p);

   //  t1 = _mm256_slli_epi32(t3, 1);

   //  t1 = _mm256_srli_epi32(t1, 1);
  

  t2 = _mm256_srli_epi32(t3, MRSN_PRIME_POW); //high part
  
  t3 = _mm256_add_epi32(t1, t2);


  return t3;
  
  //Overflow? Only if a+b = 0xFFFFFFFF = 2*p+1.
  //But a<=p, b<=p, overflow here never happens
  // res = _mm256_and_si256(t3, p);
  
  // return res;
  
}



static inline __attribute__((always_inline))
__m256i sub(__m256i a, __m256i b)
{

  __m256i t1, t2, t3;

  const __m256i p = _mm256_set1_epi32(MRSN_PRIME);
    
  __m256i minus_b = _mm256_xor_si256(b, p);
  
  t3 = _mm256_add_epi32(a, minus_b);
  
  t1 = _mm256_and_si256(t3, p);

  t2 = _mm256_srli_epi32(t3, MRSN_PRIME_POW);
  
  t3  = _mm256_add_epi32(t1, t2);

  return t3;
  // res = _mm256_and_si256(t3, p);
  
  //return res;
  
}




static inline __attribute__((always_inline))
__m256i mul(__m256i a, __m256i b)
{
  // const __m256i p = _mm256_set1_epi32(MRSN_PRIME);
  
  // const __m256i lo_mask = _mm256_setr_epi32(-1,  0, -1,  0, -1,  0, -1,  0); 

    
  __m256i prod_even = _mm256_mul_epu32(a, b);  //[_ a6 _ a4 _ a2 _ a0] .* [_ b6 _ b4 _ b2 _ b0] = [ab6  ab4  ab2  ab0]
                                               // 8 x 32                   8 x 32                  4 x 64 

  __m256i a_odd = _mm256_srli_epi64(a, 32);    //[0 a7 0 a5 0 a3 0 a1]
  __m256i b_odd = _mm256_srli_epi64(b, 32);    //[0 b7 0 b5 0 b3 0 b1]
  
  __m256i prod_odd = _mm256_mul_epu32(a_odd, b_odd); //[0 a7 0 a5 0 a3 0 a1] .* [0 b7 0 b5 0 b3 0 b1] = [ab7  ab5  ab3  ab1]
                                                     //                                                  4 x 64

  //ab = ab_hi * 2^31 + ab_lo, each with 31 bits

  __m256i hi_even = _mm256_srli_epi64(prod_even, MRSN_PRIME_POW); //[0 ab6_hi 0 ab4_hi 0 ab2_hi 0 ab0_hi]

  __m256i temp    = _mm256_srli_epi64(prod_odd, MRSN_PRIME_POW);  //[0 ab7_hi 0 ab5_hi 0 ab3_hi 0 ab1_hi]
  
  __m256i hi_odd  = _mm256_slli_epi64(temp, 32);                      //[ab7_hi 0 ab5_hi 0 ab3_hi 0 ab1_hi 0]
  
  __m256i hi      = _mm256_or_si256(hi_even, hi_odd);
  

  __m256i lo_odd  = _mm256_slli_epi64(prod_odd, 33);  //[ab7_lo 0  ab5_lo 0  ab3_lo  0 ab1_lo 0]

  lo_odd  =   _mm256_srli_epi64(lo_odd, 1);
    
  __m256i temp2   = _mm256_slli_epi64(prod_even, 33); //[ab6_lo 0  ab4_lo 0  ab2_lo  0 ab0_lo 0]

  __m256i lo_even = _mm256_srli_epi64(temp2, 33); //[0 ab6_lo 0  ab4_lo 0  ab2_lo  0 ab0_lo]

  __m256i lo      = _mm256_or_si256(lo_even, lo_odd); 
    
  return  add(lo, hi);
  
}




static inline __attribute__((always_inline))
__m256i rot15(__m256i a)
{
  //   return ((x << 15) | (x >> 16)) & 0x7FFFFFFF;

  __m256i t1, t2, t3, res;

  const __m256i p = _mm256_set1_epi32(MRSN_PRIME);
   
  t1 = _mm256_slli_epi32(a, 15);

  t2 = _mm256_srli_epi32(a, 16);

  t3 = _mm256_or_si256(t1, t2);

  res =  _mm256_and_si256(t3, p);

  return res;
}



static inline __attribute__((always_inline))
__m256i rot(__m256i a, int k)
{
  //   return ((x << k) | (x >> 31-k)) & 0x7FFFFFFF;

  __m256i t1, t2, t3, res;

  int km =  MRSN_PRIME_POW - k;

  const __m256i p = _mm256_set1_epi32(MRSN_PRIME);
   
  t1 = _mm256_slli_epi32(a, k);

  t2 = _mm256_srli_epi32(a, km);

  t3 = _mm256_or_si256(t1, t2);

  res =  _mm256_and_si256(t3, p);

  return res;
}

static inline __attribute__((always_inline))
void bf1(__m256i * a_re, __m256i * a_im, __m256i * b_re, __m256i * b_im)
{
  /*
    a = a + b
    b = a - b
  */

  __m256i t1, t2;
  
  t1    = add(*a_re, *b_re);
  *b_re = sub(*a_re, *b_re);
  *a_re = t1;

  t2    = add(*a_im, *b_im);
  *b_im = sub(*a_im, *b_im);
  *a_im = t2;

}

static inline __attribute__((always_inline))
void twister(__m256i Ar,
	     __m256i Ai,
	     __m256i Br,
	     __m256i Bi,
	     __m256i c,
	     __m256i s,
	     int32_t * S0,
	     int32_t * S1,
	     int32_t * S2,
	     int32_t * S3)
{
  __m256i t0, t1, t2, AipBr, AimBr, ArpBi, ArmBi;

 
  AipBr = add(Ai, Br);					
  AimBr = sub(Ai, Br);					
  ArpBi = add(Ar, Bi);					
  ArmBi = sub(Ar, Bi);					
  t0 = mul(ArmBi, c);				
  t1 = mul(AipBr, s);				
  t2 = sub(t0, t1);
  _mm256_store_si256((__m256i *)S0, t2);					
  t0 = mul(AipBr, c);				
  t1 = mul(ArmBi, s);				
  t2 = add(t0, t1);
   _mm256_store_si256((__m256i *)S1, t2);					
  t0 = mul(ArpBi, c);				
  t1 = mul(AimBr, s);				
  t2 = add(t0, t1);
  _mm256_store_si256((__m256i *)S2, t2);					
  t0 = mul(AimBr, c);				
  t1 = mul(ArpBi, s);				
  t2 = sub(t0, t1);
  _mm256_store_si256((__m256i *)S3, t2);

  return;
};


static inline __attribute__((always_inline))
void twister2(__m256i Ar,
	     __m256i Ai,
	     __m256i Br,
	     __m256i Bi,
	     __m256i c,
	     __m256i s,
	     int32_t * S0,
	     int32_t * S1,
	     int32_t * S2,
	     int32_t * S3)
{
  __m256i t0, t1, t2;// AipBr, AimBr, ArpBi, ArmBi;

 
  //AipBr = add(Ai, Br);					
  //AimBr = sub(Ai, Br);
  add_sub(&Ai, &Br);
  
  //ArpBi = add(Ar, Bi);					
  //ArmBi = sub(Ar, Bi);
  add_sub(&Ar, &Bi);
  
  //t0 = mul(ArmBi, c);				
  //t1 = mul(AipBr, s);
  t0 = mul(Bi, c);				
  t1 = mul(Ai, s);
  
  t2 = sub(t0, t1);
  _mm256_store_si256((__m256i *)S0, t2);
  
  //t0 = mul(AipBr, c);				
  //t1 = mul(ArmBi, s);
  t0 = mul(Ai, c);				
  t1 = mul(Bi, s);
    
  t2 = add(t0, t1);
   _mm256_store_si256((__m256i *)S1, t2);

   
   //t0 = mul(ArpBi, c);				
   //t1 = mul(AimBr, s);
   t0 = mul(Ar, c);				
   t1 = mul(Br, s);
   t2 = add(t0, t1);
  _mm256_store_si256((__m256i *)S2, t2);
  
  //t0 = mul(AimBr, c);				
  //t1 = mul(ArpBi, s);
  t0 = mul(Br, c);				
  t1 = mul(Ar, s);
  t2 = sub(t0, t1);
  _mm256_store_si256((__m256i *)S3, t2);

  return;
};


/*
static inline __attribute__((always_inline))
void twister3(__m256i Ar,
	     __m256i Ai,
	     __m256i Br,
	     __m256i Bi,
	     __m256i c,
	     __m256i s,
	     int32_t * S0,
	     int32_t * S1,
	     int32_t * S2,
	     int32_t * S3)
{
  __m256i cc, t0, t1, t2, AipBr, AimBr, ArpBi, ArmBi;

  AipBr = add(Ai, Br);					
  AimBr = sub(Ai, Br);					


  ArpBi = add(Ar, Bi);					
  ArmBi = sub(Ar, Bi);					



  
  t0 = mul(ArmBi, c);				
  t1 = mul(AipBr, s);				
  t2 = sub(t0, t1);
  _mm256_store_si256((__m256i *)S0, t2);					
  t0 = mul(AipBr, c);				
  t1 = mul(ArmBi, s);				
  t2 = add(t0, t1);
   _mm256_store_si256((__m256i *)S1, t2);					
  t0 = mul(ArpBi, c);				
  t1 = mul(AimBr, s);				
  t2 = add(t0, t1);
  _mm256_store_si256((__m256i *)S2, t2);					
  t0 = mul(AimBr, c);				
  t1 = mul(ArpBi, s);				
  t2 = sub(t0, t1);
  _mm256_store_si256((__m256i *)S3, t2);


  cc = c;

  add_sub(&c, &s);

  
  
  //AipBr = add(Ai, Br);					
  //AimBr = sub(Ai, Br);
  add_sub(&Ai, &Br);
  
  //ArpBi = add(Ar, Bi);					
  //ArmBi = sub(Ar, Bi);
  add_sub(&Ar, &Bi);
  
  //t0 = mul(ArmBi, c);				
  //t1 = mul(AipBr, s);
  t0 = mul(Bi, c);				
  t1 = mul(Ai, s);
  
  t2 = sub(t0, t1);
  _mm256_store_si256((__m256i *)S0, t2);
  
  //t0 = mul(AipBr, c);				
  //t1 = mul(ArmBi, s);
  t0 = mul(Ai, c);				
  t1 = mul(Bi, s);
    
  t2 = add(t0, t1);
   _mm256_store_si256((__m256i *)S1, t2);

   
   //t0 = mul(ArpBi, c);				
   //t1 = mul(AimBr, s);
   t0 = mul(Ar, c);				
   t1 = mul(Br, s);
   t2 = add(t0, t1);
  _mm256_store_si256((__m256i *)S2, t2);
  
  //t0 = mul(AimBr, c);				
  //t1 = mul(ArpBi, s);
  t0 = mul(Br, c);				
  t1 = mul(Ar, s);
  t2 = sub(t0, t1);
  _mm256_store_si256((__m256i *)S3, t2);

  return;
};


*/

#endif

