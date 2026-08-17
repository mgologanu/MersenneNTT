
extern const uint64_t rad8_d32[];
extern const uint64_t rad8_d64[];
extern const uint64_t rad8_d128[];
extern const uint64_t d16[];
extern const uint64_t d32[];
extern const uint64_t d64[];
extern const uint64_t d128[];
extern const uint64_t d256[];
extern const uint64_t d512[];
extern const uint64_t d1024[];
extern const uint64_t d2048[];
extern const uint64_t d4096[];
extern const uint64_t d8192[];
extern const uint64_t d16384[];
extern const uint64_t d32768[];
extern const uint64_t d65536[];
extern const uint64_t d131072[];
extern const uint64_t d262144[];
extern const uint64_t d524288[];
extern const uint64_t d1048576[];


extern void cspass_sr(uint64_t *a,  const uint64_t *w,  size_t n);
extern void cspass_srl3(uint64_t *a,  const uint64_t *w,  const uint64_t *w2, size_t n);

extern void cspass_rad4(uint64_t *a,  const uint64_t *w,  size_t n);
extern void cspass_rad8(uint64_t *a,  const uint64_t *w, size_t n);

extern void cs16_vert(uint64_t *a);
extern void cs8_vert(uint64_t *a);
extern void c4_vert(uint64_t *a);

extern void ts16_V4(uint64_t *a);




extern void uspass_sr(uint64_t *a,  const uint64_t *w,  size_t n);
extern void uspass_srl3(uint64_t *a,  const uint64_t *w,  const uint64_t *w2, size_t n);

extern void uspass_rad4(uint64_t *a,  const uint64_t *w,  size_t n);
extern void uspass_rad8(uint64_t *a,  const uint64_t *w, size_t n);

extern void us16_vert(uint64_t *a);
extern void us8_vert(uint64_t *a);
extern void u4_vert(uint64_t *a);

extern void uts16_V4(uint64_t *a);

extern void rpass_sr(uint64_t *a,  const uint64_t *w,  size_t n);
extern void rpass_srl3(uint64_t *a,  const uint64_t *w,  const uint64_t *w2, size_t n);
extern void rpass_srl3_half(uint64_t *a,  const uint64_t *w,  const uint64_t *w2, size_t n);


extern void vpass_sr(uint64_t *a,  const uint64_t *w,  size_t n);
extern void vpass_srl3(uint64_t *a,  const uint64_t *w,  const uint64_t *w2, size_t n);



static inline __attribute__((always_inline))
void red1(__m256i *T)
{
  const __m256i p = BROADCAST(PRIME);
  
  const __m256i p_64 = BROADCAST(PRIME_64);

  __m256i t1, t2, t3;

  t3  = _mm256_add_epi64(*T,  p_64);			
  t1  = _mm256_and_si256(t3, p);				
  t2  = _mm256_srli_epi64(t3, P_POW);			
  *T  =  _mm256_add_epi64( t1,  t2);
  return;
};

static inline __attribute__((always_inline))
__m256i mul1(__m256i a, __m256i b)
{
  
  const __m256i p = BROADCAST(PRIME);

  const __m256i p_64 = BROADCAST(PRIME_64);

  __m256i t1, t2, t3;
    
  t3  = _mm256_add_epi64(a, p_64);				
  t1  = _mm256_and_si256(t3, p);				
  t2  = _mm256_srli_epi64(t3, P_POW);			
  t3  = _mm256_add_epi64(t1, t2);			
  t1  = _mm256_mul_epu32(t3, b);				
  t2  = _mm256_and_si256(t1, p);				
  t3  = _mm256_srli_epi64(t1, P_POW);			
  return _mm256_add_epi64(t2, t3);				
};

static inline __attribute__((always_inline))
__m256i mul0(__m256i a, __m256i b)
{
  
  const __m256i p = BROADCAST(PRIME);

   __m256i t1, t2, t3;

  t1  = _mm256_mul_epu32(a, b);					
  t2  = _mm256_and_si256(t1, p);				
  t3  = _mm256_srli_epi64(t1, P_POW);			
  return  _mm256_add_epi64(t2, t3);			 
}


static inline __attribute__((always_inline))
void twisted(__m256i Ar,
	     __m256i Ai,
	     __m256i Br,
	     __m256i Bi,
	     __m256i c,
	     __m256i s,
	     uint64_t * S0,
	     uint64_t * S1,
	     uint64_t * S2,
	     uint64_t * S3)
{
  __m256i t0, t1, t2, AipBr, AimBr, ArpBi, ArmBi;

 
  AipBr = _mm256_add_epi64(Ai, Br);					
  AimBr = _mm256_sub_epi64(Ai, Br);					
  ArpBi = _mm256_add_epi64(Ar, Bi);					
  ArmBi = _mm256_sub_epi64(Ar, Bi);					
  t0 = mul1(ArmBi, c);				
  t1 = mul1(AipBr, s);				
  t2 = _mm256_sub_epi64(t0, t1);
  red1(&t2);
  _mm256_store_si256((__m256i *)S0, t2);					
  t0 = mul1(AipBr, c);				
  t1 = mul1(ArmBi, s);				
  t2 = _mm256_add_epi64(t0, t1);
  red1(&t2);
   _mm256_store_si256((__m256i *)S1, t2);					
  t0 = mul1(ArpBi, c);				
  t1 = mul1(AimBr, s);				
  t2 = _mm256_add_epi64(t0, t1);
  red1(&t2);
   _mm256_store_si256((__m256i *)S2, t2);					
  t0 = mul1(AimBr, c);				
  t1 = mul1(ArpBi, s);				
  t2 = _mm256_sub_epi64(t0, t1);
  red1(&t2);
  _mm256_store_si256((__m256i *)S3, t2);

  return;
};


static inline __attribute__((always_inline))
void transpose4_from_mem(__m256i * x0, __m256i * x1, __m256i * x2, __m256i * x3, uint64_t *a, size_t offset, size_t stride)
{
  __m256i t0, t1, t2,t3;
  
  t0 = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_load_si128((__m128i *) &a[offset+0*stride+0])),_mm_load_si128((__m128i *) &a[offset+2*stride+0]), 1);
  t1 = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_load_si128((__m128i *) &a[offset+1*stride+0])), _mm_load_si128((__m128i *) &a[offset+3*stride+0]), 1); 
  t2 = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_load_si128((__m128i *) &a[offset+0*stride+2])), _mm_load_si128((__m128i *) &a[offset+2*stride+2]), 1); 
  t3 = _mm256_insertf128_si256(_mm256_castsi128_si256(_mm_load_si128((__m128i *) &a[offset+1*stride+2])), _mm_load_si128((__m128i *) &a[offset+3*stride+2]), 1); 
  *x0 = _mm256_unpacklo_epi64(t0,t1);
  *x1 = _mm256_unpackhi_epi64(t0,t1);
  *x2 = _mm256_unpacklo_epi64(t2,t3);
  *x3 = _mm256_unpackhi_epi64(t2,t3);
}


static inline __attribute__((always_inline))
void transpose4(__m256i * x0, __m256i * x1, __m256i * x2, __m256i * x3)
{
  __m256i s0, s1, s2,s3;
  
  s0 =  (__m256i ) _mm256_shuffle_pd((__m256d) (*x0) ,(__m256d)(*x1),0x0) ;					   
  s1 =  (__m256i ) _mm256_shuffle_pd((__m256d) (*x0) ,(__m256d)(*x1),0xF) ;                                       
  s2 =  (__m256i ) _mm256_shuffle_pd((__m256d) (*x2) ,(__m256d)(*x3),0x0) ;                                       
  s3 =  (__m256i ) _mm256_shuffle_pd((__m256d) (*x2) ,(__m256d)(*x3),0xF) ;					    
  *x0 = _mm256_permute2f128_si256(s0,s2,0x20) ;                                       
  *x2 = _mm256_permute2f128_si256(s0,s2,0x31) ;                                       
  *x1 = _mm256_permute2f128_si256(s1,s3,0x20) ;                                       
  *x3 = _mm256_permute2f128_si256(s1,s3,0x31) ;                                       
};





static inline __attribute__((always_inline))
void c4_reg_to_mem(__m256i x0, __m256i y0, __m256i x1, __m256i y1, __m256i x2, __m256i y2, __m256i x3, __m256i y3,
		   uint64_t *X0, uint64_t * Y0, uint64_t * X1, uint64_t * Y1, uint64_t * X2, uint64_t * Y2, uint64_t * X3, uint64_t * Y3)
{

  __m256i t0, t1, t2, t3;
  __m256i AipBr, AimBr, ArpBi, ArmBi;
  __m256i Ar, Br, Ai, Bi, t02, t13, s02, s13;
  
  t02 = _mm256_add_epi64(x0, x2);          /* x 0+2 */				
  Ar  = _mm256_sub_epi64(x0, x2);          /* x 0-2 */				
  t13 = _mm256_add_epi64(x1, x3);          /* x 1+3 */				
  Br  = _mm256_sub_epi64(x1, x3);          /* x 1-3 */				
  s02 = _mm256_add_epi64(y0, y2);          /* y 0+2 */				
  Ai  = _mm256_sub_epi64(y0, y2);          /* y 0-2 */				
  s13 = _mm256_add_epi64(y1, y3);          /* y 1+3 */				
  Bi  = _mm256_sub_epi64(y1, y3);          /* y 1-3 */				
  t0  = _mm256_add_epi64(t02, t13);         /* x 0+2+1+3 */			
  t1  = _mm256_add_epi64(s02, s13);         /* y 0+2+1+3 */			
  t2  = _mm256_sub_epi64(t02, t13);         /* x 0+2-1-3 */			
  t3  = _mm256_sub_epi64(s02, s13);         /* y 0+2-1-3 */			
  AipBr = _mm256_add_epi64(Ai, Br);        /* y 0-2 + x 1-3 */			
  AimBr = _mm256_sub_epi64(Ai, Br);        /* y 0-2 - x 1-3 */			
  ArpBi = _mm256_add_epi64(Ar, Bi);        /* x 0-2 + y 1-3 */			
  ArmBi = _mm256_sub_epi64(Ar, Bi);        /* x 0-2 - y 1-3 */

  red1(&t0);
  _mm256_store_si256((__m256i *)X0, t0);
  red1(&t1);
  _mm256_store_si256((__m256i *)Y0, t1);
  red1(&t2);
  _mm256_store_si256((__m256i *)X1, t2);
  red1(&t3);
  _mm256_store_si256((__m256i *)Y1, t3);
  red1(&ArmBi);
  _mm256_store_si256((__m256i *)X2, ArmBi);
  red1(&AipBr);
  _mm256_store_si256((__m256i *)Y2, AipBr);
  red1(&ArpBi);
  _mm256_store_si256((__m256i *)X3, ArpBi);
  red1(&AimBr);
  _mm256_store_si256((__m256i *)Y3, AimBr);							
}




static inline __attribute__((always_inline))
void c2_reg_to_mem(__m256i x0, __m256i y0, __m256i x1, __m256i y1,
		   uint64_t *X0, uint64_t * Y0, uint64_t * X1, uint64_t * Y1)
{
  __m256i t0, t1, t2, t3;

  t0 = _mm256_add_epi64(x0, x1);				     
  t2 = _mm256_sub_epi64(x0, x1);				     
  t1 = _mm256_add_epi64(y0, y1);				     
  t3 = _mm256_sub_epi64(y0, y1);				     
  
  red1(&t0);
  red1(&t1);
  red1(&t2);
  red1(&t3);
    
  _mm256_store_si256((__m256i *)X0, t0);				     
  _mm256_store_si256((__m256i *)Y0, t1);				     
  _mm256_store_si256((__m256i *)X1, t2);				     
  _mm256_store_si256((__m256i *)Y1, t3);				     
} 




/*

  Inverse of TWISTED: from (z0, z1) to (A, B):


  z0 = (A+Bj)*(c+js) = Ac - Bs + j(As + Bc)
  z1 = (A-Bj)*(c-js)

  where A and B are complex:
  
  A = Ar + j Ai
  B = Br + j Bi

  and z0 and z1 are given by:


  z0 = S0 + j S1
  z1 = S2 + j S3


  Multiply by c+/-js:

  z0*(c-js) = A+B*j
  z1*(c+js) = A-B*j


  (z0+z1)*c + j*(z1-z0)*s = 2A (forget the 2!) = Ar + jAi
  (z0-z1)*c - j*(z0+z1)*s = 2*B*j (forget the 2) = (Br + jBi)*j = -Bi + jBr

  Multiply last eq by j:

  j*(z0-z1)*c - j^2 * (z0+z1)*s = -j*Bi + j^2 Br

  j*(z0-z1)*c  +  (z0+z1)*s = -j*Bi + - Br
  
  and again by -1:

  -(z0+z1)*s + j*(z1-z0)*c = Br + j Bi = B!!

  Final eqs:
  
  (z0+z1)*c + j*(z1-z0)*s = A
 -(z0+z1)*s + j*(z1-z0)*c = B
  
  
  

*/


static inline __attribute__((always_inline))
void untwisted(__m256i *Ar,
	       __m256i *Ai,
	       __m256i *Br,
	       __m256i *Bi,
	       __m256i c,
	       __m256i s,
	       uint64_t * S0,
	       uint64_t * S1,
	       uint64_t * S2,
	       uint64_t * S3)
{
  __m256i t0, t1, z0r, z0i, z1r, z1i, t0r, t1r, t0i, mt1i;
  

  z0r = _mm256_load_si256((__m256i *) S0);						
  z0i = _mm256_load_si256((__m256i *) S1);						
  z1r = _mm256_load_si256((__m256i *) S2);						
  z1i = _mm256_load_si256((__m256i *) S3);						
  t0r = _mm256_add_epi64(z0r, z1r);						
  t0i = _mm256_add_epi64(z0i, z1i);						
  t1r = _mm256_sub_epi64(z1r, z0r);						
  mt1i= _mm256_sub_epi64(z0i, z1i);						
  t0  = mul1(t0r,c);						
  t1  = mul1(mt1i,s);						
  *Ar  = _mm256_add_epi64(t0,t1);						
  t0  = mul1(t0i,c);						
  t1  = mul1(t1r,s);						
  *Ai  = _mm256_add_epi64(t0, t1);						
  t0  = mul1(mt1i,c);						
  t1  = mul1(t0r,s);						
  *Br  = _mm256_sub_epi64(t0,t1);						
  t0  = mul1(t1r,c);						
  t1  = mul1(t0i,s);						
  *Bi  = _mm256_sub_epi64(t0,t1);						
};




static inline __attribute__((always_inline))
void u4_reg_to_mem(
		   uint64_t *X0, uint64_t *Y0,
		   uint64_t *X1, uint64_t *Y1,
		   uint64_t *X2, uint64_t *Y2,
		   uint64_t *X3, uint64_t *Y3,
		   __m256i *x0, __m256i *y0,
		   __m256i *x1, __m256i *y1,
		   __m256i *x2, __m256i *y2,
		   __m256i *x3, __m256i *y3)

{

  __m256i Ar, Br, Ai, Bi,  AipBr, AimBr, ArpBi, ArmBi, t0, t1, t2, t3;
  __m256i t02, t13, s02, s13;
 
  
  t0 = _mm256_load_si256((__m256i *)X0);
  t1 = _mm256_load_si256((__m256i *)Y0);							
  t2 = _mm256_load_si256((__m256i *)X1);							
  t3 = _mm256_load_si256((__m256i *)Y1);							
  ArmBi = _mm256_load_si256((__m256i *)X2);							
  AipBr = _mm256_load_si256((__m256i *)Y2);							
  ArpBi = _mm256_load_si256((__m256i *)X3);							
  AimBr = _mm256_load_si256((__m256i *)Y3);							
  Ai  = _mm256_add_epi64(AipBr, AimBr);	/* y 0-2 */				
  Br  = _mm256_sub_epi64(AipBr, AimBr);	/* x 1-3 */				
  Ar  = _mm256_add_epi64(ArpBi, ArmBi);	/* x 0-2 */				
  Bi  = _mm256_sub_epi64(ArpBi, ArmBi);	/* y 1-3 */				
  t02 = _mm256_add_epi64(t0, t2);            /* x 0+2 */				
  t13 = _mm256_sub_epi64(t0, t2);            /* x 1+3 */				
  s02 = _mm256_add_epi64(t1, t3);            /* y 0+2 */				
  s13 = _mm256_sub_epi64(t1, t3);            /* y 1+3 */				
  *x0  = _mm256_add_epi64(t02, Ar);							
  *x2  = _mm256_sub_epi64(t02, Ar);							
  *x1  = _mm256_add_epi64(t13, Br);							
  *x3  = _mm256_sub_epi64(t13, Br);							
  *y0  = _mm256_add_epi64(s02, Ai);							
  *y2  = _mm256_sub_epi64(s02, Ai);							
  *y1  = _mm256_add_epi64(s13, Bi);							
  *y3  = _mm256_sub_epi64(s13, Bi);							
}



static inline __attribute__((always_inline))
void u2_reg_to_mem(
		   uint64_t *X0, uint64_t *Y0,
		   uint64_t *X1, uint64_t *Y1,
		   __m256i *x0, __m256i *y0,
		   __m256i *x1, __m256i *y1)
{
  __m256i t0, t1, t2, t3;

  t0 = _mm256_load_si256((__m256i *)X0);				      
  t1 = _mm256_load_si256((__m256i *)Y0);				      
  t2 = _mm256_load_si256((__m256i *)X1);				      
  t3 = _mm256_load_si256((__m256i *)Y1);				      
  *x0 = _mm256_add_epi64(t0, t2);				      
  *x1 = _mm256_sub_epi64(t0, t2);				      
  *y0 = _mm256_add_epi64(t1, t3);				      
  *y1 = _mm256_sub_epi64(t1, t3);				      
}

static inline __attribute__((always_inline))
void untwisted_reg(__m256i *Ar,
		   __m256i *Ai,
		   __m256i *Br,
		   __m256i *Bi,
		   __m256i c,
		   __m256i s,
		   __m256i z0r,
		   __m256i z0i,
		   __m256i z1r,
		   __m256i z1i)
{
  __m256i t0, t1, t0r, t1r, t0i, mt1i;
 
  t0r = _mm256_add_epi64(z0r, z1r);						
  t0i = _mm256_add_epi64(z0i, z1i);						
  t1r = _mm256_sub_epi64(z1r, z0r);						
  mt1i= _mm256_sub_epi64(z0i, z1i);						
  t0  = mul1(t0r,c);						
  t1  = mul1(mt1i,s);						
  *Ar  = _mm256_add_epi64(t0,t1);						
  t0  = mul1(t0i,c);						
  t1  = mul1(t1r,s);						
  *Ai  = _mm256_add_epi64(t0, t1);						
  t0  = mul1(mt1i,c);						
  t1  = mul1(t0r,s);						
  *Br  = _mm256_sub_epi64(t0,t1);						
  t0  = mul1(t1r,c);						
  t1  = mul1(t0i,s);						
  *Bi  = _mm256_sub_epi64(t0,t1);						
}




static inline __attribute__((always_inline))
void real_twisted(__m256i Ar,
		  __m256i Br,
		  __m256i c,
		  __m256i s,
		  uint64_t * S0,
		  uint64_t * S1)
{
  __m256i t0, t1, t2;

  t0 = mul1(Ar, c);				
  t1 = mul1(Br, s);				
  t2 = _mm256_sub_epi64(t0, t1);
  red1(&t2);
  _mm256_store_si256((__m256i *)S0, t2);
  
  t0 = mul1(Br, c);				
  t1 = mul1(Ar, s);				
  t2 = _mm256_add_epi64(t0, t1);
  red1(&t2);
  _mm256_store_si256((__m256i *)S1, t2);

  
}



static inline __attribute__((always_inline))
void real_untwisted(__m256i * Ar,
		    __m256i * Br,
		    __m256i c,
		    __m256i s,
		    uint64_t * S0,
		    uint64_t * S1)
{
  __m256i t0, t1, t2, t3;

  t0 = _mm256_load_si256((__m256i *)S0);				
  t1 = _mm256_load_si256((__m256i *)S1);				
  t2 = mul1(t0,c);				
  t3 = mul1(t1,s);				
  *Ar = _mm256_add_epi64(t2,t3);
  red1(Ar);
  t2 = mul1(t0,s);				
  t3 = mul1(t1,c);				
  *Br = _mm256_sub_epi64(t3,t2);				
  red1(Br);
}



static inline __attribute__((always_inline))
void red1_scalar(uint64_t *T)
{
  const uint64_t p = PRIME;
  
  const uint64_t p_64 = PRIME_64;

  uint64_t t1, t2, t3;

  t3  = *T +  p_64;			
  t1  = t3 & p;				
  t2  = t3 >> P_POW;			
  *T  =  t1 + t2;
  return;
};

