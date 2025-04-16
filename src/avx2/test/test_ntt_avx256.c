#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "minunit.h"

#include<immintrin.h>
#include <stdalign.h>

/* typedef int64_t integer; */
/* typedef __m512i integer_packed; */

#include "ntt.h"

#include "pre8_V4.h"

#include "i.h"

integer_packed mul_mod_mersenne_avx256_64(integer_packed va, integer_packed vb, integer_packed vp);
integer_packed mul_mod_mersenne_avx256_64_2(integer_packed va, integer_packed vb, integer_packed vp);
integer_packed red_mersenne_avx256_64(integer_packed vT, integer_packed vp);
integer_packed red_mersenne_avx256_64_2(integer_packed vT, integer_packed vp);
integer_packed red_mersenne_avx256_64_2_final(integer_packed vT, integer_packed vp);


MU_TEST(test_mul_mod_avx256) {

  // int V = 8;

  uint64_t a;
  uint64_t b;
  
  uint64_t T;
  
  uint64_t r_expected;

  uint64_t p = PRIME;
  
  uint64_t *va = (uint64_t *) aligned_alloc(64, V*sizeof(uint64_t));
  uint64_t *vb = (uint64_t *) aligned_alloc(64, V*sizeof(uint64_t));
  uint64_t *vr = (uint64_t *) aligned_alloc(64, V*sizeof(uint64_t));
  uint64_t *vp = (uint64_t *) aligned_alloc(64, V*sizeof(uint64_t));
 

  integer_packed pa;
  integer_packed pb;
  integer_packed pp;
  integer_packed pr;

  integer_packed t1,t2,t3;

  for (int i = 0; i<V; i++)
    {
      vp[i] = p;
    }
    
  pp = _mm256_load_si256((__m256i *) vp);


  
  a = 2*(p-1);
  b = 2*(p-1)+3;

  T = a*b;

  // printf("%lu %lu %lu\n", a,b,T);
  
  for (int i = 0; i<V; i++)
    {
      va[i] = a;
      vb[i] = b;
    }
  
  pa = _mm256_load_si256((__m256i *)va);
  pb = _mm256_load_si256((__m256i *)vb);

  pr = _mm256_mul_epu32(pa, pb);

  _mm256_store_si256((__m256i *) vr, pr);

   for (int i = 0; i<V; i++)
    {
      //      printf("%lu ", vr[i]);
      mu_assert(vr[i] == T, "Mul gresit");
    }
   //  printf("\n");

   pr = _mm256_sub_epi64(pa, pb);

   
   _mm256_store_si256((__m256i *) vr, pr);

   r_expected = 18446744073709551613UL;
   
   for (int i = 0; i<V; i++)
    {
      //      printf("%lu ", vr[i]);
      mu_assert(vr[i] == r_expected, "Mul gresit");
    }
   //  printf("\n");



  
   a = 18446744073709551613UL; //2^64-3
   b = 3;

   T = a*b;

   //   printf("%lu %lu %lu\n", a,b,T);

   r_expected = 18446744073709551607UL; //2^64-9
   
   mu_assert(T == r_expected, "Error unsigned scalar multiplication");

  
   for (int i = 0; i<V; i++)
     {
       va[i] = a;
       vb[i] = b;
     }
  
   pa = _mm256_load_si256((__m256i *)va);
   pb = _mm256_load_si256((__m256i *)vb);
   
   pr = _mm256_mul_epu32(pa, pb);
   
   _mm256_store_si256((__m256i *) vr, pr);
   
   for (int i = 0; i<V; i++)
    {
      //      printf("%lu ", vr[i]);
      mu_assert(vr[i] != T, "Vectorial multiplication is not the same, as only the smallest 32 bit are multiplied");
    }
   // printf("\n");


   /*
     The value of a can both increase due to lazy evaluation and be wround up with respect to 2^64 due to
      difference ( if x-y < 0 -> 2^64 + (x-y)).

      We need to reduce a below 2^32-1 so that multiplication of
      lowest 32 bits can be applied; it is enough to apply 2x reduction
     method to do reduction also for "negative" values, by adding a large multiple of p, here 1024*p
     and then do a single reduction step to arrive in [0..2p] < 2^32-1

   */
  uint64_t p_x_1024 = 2199023254528UL;

  t2  = _mm256_set1_epi64x(p_x_1024);

  t3 = _mm256_add_epi64(pa, t2);

  t1 = _mm256_and_si256(t3, pp);
 
  t2 = _mm256_srli_epi64(t3, 31);
  
  t3 = _mm256_add_epi64(t1, t2);

  
   
   
   pr = mul_mod_mersenne_avx256_64(t3, pb, pp);

   _mm256_store_si256((__m256i *) vr, pr);

   r_expected = 2147483638UL; //p-9;
   // printf("\n");
   for (int i = 0; i<V; i++)
    {
      //      printf("%lu ", vr[i]);
      mu_assert(vr[i] == r_expected, "Modular multiplication with avx256 and negative a did not work");
    }
   //  printf("\n");


  a = 2*(p+1)-1;
  b = 2*(p+1)-1;
  for (int i = 0; i<V; i++)
    {
      va[i] = a;
      vb[i] = b;
    }
  
  pa = _mm256_load_si256((__m256i *)va);
  pb = _mm256_load_si256((__m256i *)vb);
  
  r_expected = p+1;
  pr = mul_mod_mersenne_avx256_64_2(pa, pb, pp);
  _mm256_store_si256((__m256i *) vr, pr);
  for (int i = 0; i<V; i++)
    {
      mu_assert(vr[i] == r_expected, "Mul gresit");
      //     printf("%ld ", vr[i]);
    }
  //  printf("\n");

   
  free(va);
  free(vb);
  free(vr);
  free(vp);

}

MU_TEST(test_reduction_avx256) {

  // int V = 8;

  uint64_t a;
  uint64_t b;

  uint64_t p = PRIME;
  uint64_t T;
  
  uint64_t r_expected;

  
  uint64_t *vT = (uint64_t *) aligned_alloc(64, V*sizeof(uint64_t));
  uint64_t *vr = (uint64_t *) aligned_alloc(64, V*sizeof(uint64_t));
  uint64_t *vp = (uint64_t *) aligned_alloc(64, V*sizeof(uint64_t));


  integer_packed pT;
  integer_packed pp;
  integer_packed pr;
 

  for (int i = 0; i<V; i++)
    {
      vp[i] = p;
 
    }
    
  pp        = _mm256_load_si256((__m256i *) vp);

  int64_t pppp = p;
  
  pp = _mm256_set1_epi64x(pppp);

  a = 2*(p+1)-1;
  b = 2*(p+1)-1;
  T = a*b;

  
  for (int i = 0; i<V; i++)
    {
      vT[i] = T;
    }
  pT = _mm256_load_si256((__m256i *)vT);


  
  r_expected = 8589934589;
  pr = red_mersenne_avx256_64(pT, pp);
  
  _mm256_store_si256((__m256i *) vr, pr);
  for (int i = 0; i<V; i++)
    {
      //   printf("%ld ", vr[i]);
      mu_assert(vr[i] == r_expected, "Error in avx256 modular reduction");
     }
  //  printf("\n");
  
  
  
  r_expected = p+1;
  pr = red_mersenne_avx256_64_2(pT, pp);
  
  _mm256_store_si256((__m256i *) vr, pr);
  for (int i = 0; i<V; i++)
    {
      //    printf("%ld ", vr[i]);
      mu_assert(vr[i] == r_expected, "Error in avx256 modular reduction");
     }
  // printf("\n");

  
  
  r_expected = 1;
  pr = red_mersenne_avx256_64_2_final(pT, pp);
  
  _mm256_store_si256((__m256i *) vr, pr);
  for (int i = 0; i<V; i++)
    {
      //printf("%ld ", vr[i]);
	  mu_assert(vr[i] == r_expected, "Error in avx256 modular reduction");
     }
  //  printf("\n");

  T = T-1;

  
  for (int i = 0; i<V; i++)
    {
      vT[i] = T;
    }
  pT = _mm256_load_si256((__m256i *)vT);

  r_expected = p;
  pr = red_mersenne_avx256_64_2(pT, pp);
  
  _mm256_store_si256((__m256i *) vr, pr);
  for (int i = 0; i<V; i++)
    {
      //     printf("%ld ", vr[i]);
	  mu_assert(vr[i] == r_expected, "Error in avx256 modular reduction");
     }
  // printf("\n");

  
  r_expected = 0;
  pr = red_mersenne_avx256_64_2_final(pT, pp);
  
  _mm256_store_si256((__m256i *) vr, pr);
  for (int i = 0; i<V; i++)
    {
      //  printf("%ld ", vr[i]);
	 mu_assert(vr[i] == r_expected, "Error in avx256 modular reduction");
     }
  // printf("\n");
  
  
  free(vT);
  free(vr);
  free(vp);


}


MU_TEST(test_macros_avx256) {

 
  uint64_t *vr = (uint64_t *) aligned_alloc(64, V*sizeof(uint64_t));


  const integer_packed p    =  BROADCAST(PRIME);
  const integer_packed p_64 = BROADCAST(PRIME_64);
    
  integer_packed a;
  integer_packed b;

  integer_packed r;
  integer_packed _mul1_t1, _mul1_t2, _mul1_t3, retval;

  integer_packed _red1_t1, _red1_t2, _red1_t3;
  integer_packed _red3_t1, _red3_t2, _red3_t3, _red3_tk;
  

 
  a = BROADCAST(18446744073709551613UL);
  b = BROADCAST(3UL);
  
  
  uint64_t r_expected;
  
  r_expected= PRIME-9;
  
  r = MUL1(a,b);
  // _mm256_store_si256((__m256i *)vr, r);

  STORE0(vr[0], r);
  
  for (int i = 0; i<V; i++)
    {
      //      printf("%ld ", vr[i]);
       mu_assert(vr[i] == r_expected, "Error in avx256 modular reduction");
    }
  //  printf("\n");

  r_expected= PRIME-3;
  r = RED1(a);
  STORE0(vr[0], r);

  for (int i = 0; i<V; i++)
    {
      //printf("%ld ", vr[i]);
      mu_assert(vr[i] == r_expected, "Error in avx256 modular reduction");
    }
  //  printf("\n");

      
  r_expected= PRIME+3;

  a = BROADCAST(2*PRIME+3);
  r = RED1(a);
  STORE0(vr[0], r);

  for (int i = 0; i<V; i++)
    {
      //  printf("%ld ", vr[i]);
      mu_assert(vr[i] == r_expected, "Error in avx256 modular reduction");
    }
  //  printf("\n");

  r_expected= 3;
  a = BROADCAST(2*PRIME+3);
  r = RED3(a);
  STORE0(vr[0], r);

  for (int i = 0; i<V; i++)
    {
      //  printf("%ld ", vr[i]);
      mu_assert(vr[i] == r_expected, "Error in avx256 modular reduction");
    }
  //    printf("\n");


  
  r_expected= 0;
  a = BROADCAST(PRIME);
  r = RED3(a);
  STORE0(vr[0], r);

  for (int i = 0; i<V; i++)
    {
      //   printf("%ld ", vr[i]);
       mu_assert(vr[i] == r_expected, "Error in avx256 modular reduction");
    }
  //   printf("\n");
  
  
  
}
  
MU_TEST(transpose_4x4) {

  uint64_t *v = (uint64_t *) aligned_alloc(64, 4*V*sizeof(uint64_t));


  integer_packed x0;
  integer_packed x1;
  integer_packed x2;
  integer_packed x3;

  integer_packed  _t0, _t1, _t2, _t3;
  integer_packed  _s0, _s1, _s2, _s3;
  for (int i = 0; i<V; i++)
    {
      v[i] = 0;
      v[i+V] = 1;
      v[i+2*V] = 2;
      v[i+3*V] = 3;
    }
  for (int i = 0; i<4*V; i++)
    {
      //    printf("%lu ", v[i]);
    }
  //  printf("\n");
  
  TRANSPOSE4_FROM_MEM(x0,x1,x2,x3,v,0,V);
  
  STORE0(v[0], x0);
  STORE0(v[V], x1);
  STORE0(v[2*V], x2);
  STORE0(v[3*V], x3);

  for (int i = 0; i<4*V; i++)
    {
      //   printf("%lu %d   ", v[i], i-4*(i/4));
      mu_assert(v[i] == i-4*(i/4), "Error in transpose macro");
    }
  //  printf("\n");

  TRANSPOSE4_FROM_MEM(x0,x1,x2,x3,v,0,V);
  STORE0(v[0], x0);
  STORE0(v[V], x1);
  STORE0(v[2*V], x2);
  STORE0(v[3*V], x3);

  for (int i = 0; i<4*V; i++)
    {
      //  printf("%lu  ", v[i]);
      mu_assert(v[i] == i/4, "Error in transpose macro");
    }
  //   printf("\n");

  x0 = LOAD(v[0]);
  x1 = LOAD(v[V]);
  x2 = LOAD(v[2*V]);
  x3 = LOAD(v[3*V]);

  TRANSPOSE4(x0,x1,x2,x3);

  STORE0(v[0], x0);
  STORE0(v[V], x1);
  STORE0(v[2*V], x2);
  STORE0(v[3*V], x3);

  for (int i = 0; i<4*V; i++)
    {
      //  printf("%lu ", v[i]);
      mu_assert(v[i] == i-4*(i/4), "Error in transpose macro");
    }
  // printf("\n");
  
  for (int i = 0; i<4*V; i++)
    {
      v[i] = i;
    }
  /* for (int i = 0; i<4*V; i++) */
  /*   { */
  /*     printf("%lu ", v[i]); */
  /*   } */
  /* printf("\n"); */
  
  TRANSPOSE4_FROM_MEM(x0,x1,x2,x3,v,0,V);
  
  STORE0(v[0], x0);
  STORE0(v[V], x1);
  STORE0(v[2*V], x2);
  STORE0(v[3*V], x3);

  /* for (int i = 0; i<4*V; i++) */
  /*   { */
  /*     printf("%lu ", v[i]); */
  /*   } */
  /* printf("\n"); */
  
  
}


MU_TEST(rad8) {

  size_t N = 8;


  _Alignas(64) uint64_t a[] ={				\
    8,8,8,8,						\
    7,7,7,7,						\
    5,5,5,5,						\
    3,3,3,3,						\
    1,1,1,1,						\
    1,1,1,1,						\
    4,4,4,4,						\
    1,1,1,1,						\
    1,1,1,1,						\
    3,3,3,3,						\
    4,4,4,4,						\
    9,9,9,9,						\
    8,8,8,8,						\
    2,2,2,2,						\
    3,3,3,3,						\
    4,4,4,4						\
  };

  const _Alignas(32) integer w[] = {
#include "roots/p4_31_rad8_32.txt"
  };

  /* for (int i = 0; i<2*4*N; i++) */
  /*   { */
  /*     printf("%lu ", a[i]); */
  /*   } */
  /* printf("\n"); */
  
  cspass_rad8(a, w, 4);


  for (int i = 0; i<2*4*N; i++)
    {
      if (a[i]>=PRIME)
	a[i]=a[i]-PRIME;
    }
  

  /* for (int i = 0; i<2*4*N; i++) */
  /*   { */
  /*     printf("%lu ", a[i]); */
  /*   } */
  /* printf("\n"); */

 uint64_t a_octave0[] ={30,34,6,2147483645,8,6,2147483639,2147483645,2147418117,131066,65542,2147352569,2147155975,393208,327688,2147090423};
	      
  uint64_t a_octave1[] ={30,34,262144,131072,139717165,2136645213,621478439,420341104,1655850730,942453770,258735297,1846098677,717843711,473277361,843544321,1608287105};

  uint64_t a_octave2[] ={30,34,2,6,65536,458752,2147155967,196608,782874255,1023754892,628704656,423953321,32517085,1859633491,947882288,269783722};

  uint64_t a_octave3[] = {30,34,2147352575,262144,635930873,1145407406,1942732312,1906603010,1444075350,377227698,591025993,1150066190,1604789861,1619465703,1858586954,2040365251};
       
  for (int i = 0; i<2*N; i++)
    {
      //     printf("%lu ", a[4*i]);
      mu_assert(a[4*i] == a_octave0[i], "Error in cspass_rad8");
    }
  //  printf("\n");
	      
  for (int i = 0; i<2*N; i++)
    {
      //  printf("%lu ", a[4*i+1]);
      mu_assert(a[4*i+1] == a_octave1[i], "Error in cspass_rad8");
    }
  //  printf("\n");

    for (int i = 0; i<2*N; i++)
    {
      //  printf("%lu ", a[4*i+2]);
      mu_assert(a[4*i+2] == a_octave2[i], "Error in cspass_rad8");
    }
    //  printf("\n");

    for (int i = 0; i<2*N; i++)
    {
      //printf("%lu ", a[4*i+3]);
      mu_assert(a[4*i+3] == a_octave3[i], "Error in cspass_rad8");
    }
    //    printf("\n");
  

}
MU_TEST(ntt_16) {

  size_t N = 16;

  const integer_packed  p = BROADCAST(PRIME);
  
  const integer_packed p_64 = BROADCAST(PRIME_64);

  integer_packed _red3_t1, _red3_t2, _red3_t3, _red3_tk, retval;
 
 

  integer_packed pa;
   
  uint64_t *a = (uint64_t *) aligned_alloc(64, 2*N*sizeof(uint64_t));

  for (int i = 0; i<2*N; i++)
    {
      a[i] = i;
    }

  /* for (int i = 0; i<2*N; i++) */
  /*   { */
  /*     printf("%lu ", a[i]); */
  /*   } */
  /* printf("\n"); */

  _Alignas(64) uint64_t a_expected[] = {120, 524288, 1910741298, 236742333, 376, 2146959343, 236742333, 1910741298, 2147483639, 2146959359, 237790925, 1909692706, 2147483639, 524272, 1909692706, 237790925, 0, 524272, 659507844, 660556404, 2147483631, 2146959359, 1487975787, 1486927227, 2147483631, 2146959343, 1486927227, 1487975787, 0, 524288, 660556404, 659507844};

  cs16(a);

  
  /* for (int i = 0; i<2*N; i++) */
  /*   { */
  /*     printf("%lu ", a[i]); */
  /*   } */
  /* printf("\n"); */
  
  uint64_t * a1 = a;
  
  for (int i = 0; i<2*N/4; i++)
    {
      pa = _mm256_load_si256((__m256i *) a1); 
      pa = RED3(pa);
      _mm256_store_si256((__m256i *) a1, pa);
      a1 += 4;
    }

  

  for (int i = 0; i<2*N; i++)
    {
      //      printf("%lu, ", a[i]);
      mu_assert(a[i]==a_expected[i],"Error in ntt 16");
    }
  //  printf("\n");
}


MU_TEST(ntt_32) {

  size_t N = 32;

  const integer_packed  p = BROADCAST(PRIME);
  
  const integer_packed p_64 = BROADCAST(PRIME_64);

  integer_packed _red3_t1, _red3_t2, _red3_t3, _red3_tk, retval;
 
   
  integer_packed pa;
   
  uint64_t *a = (uint64_t *) aligned_alloc(64, 2*N*sizeof(uint64_t));

  for (int i = 0; i<2*N; i++)
    {
       a[i] = ((uint64_t) i) * ((uint64_t) i);
    }

  /* for (int i = 0; i<2*N; i++) */
  /*   { */
  /*     printf("%lu ", a[i]); */
  /*   } */
  /* printf("\n"); */

  _Alignas(64) uint64_t a_expected[] = {  10416, 102761536, 388816467, 1103562237, 74928, 2116024383, 865312623, 627064033, 2147483151, 2044724287, 1968384684, 850979586, 2147482127, 31455296, 1223446928, 1595915294, 1056, 98564160, 1854134033, 334228580, 2147481631, 2111831103, 523847363, 882994197, 2147481631, 2048915519, 83636718, 2006189211, 2147482655, 35650624, 1682352700, 1188990186, 828074031, 485586995, 243436932, 1581634909, 2046274808, 1039980478, 1836525048, 363432084, 1441938631, 99552763, 842886550, 2117869293, 1176730166, 595107597, 347329668, 235592063, 1957023213, 229983899, 482120574, 2111840748, 827006761, 1424575608, 42316281, 324214113, 504144155, 896157494, 948131175, 262000071, 144285095, 1335959739, 1901021065, 1392014475};
    
  cs32(a);

  
  /* for (int i = 0; i<2*N; i++) */
  /*   { */
  /*     printf("%lu ", a[i]); */
  /*   } */
  /* printf("\n"); */
  
  uint64_t * a1 = a;
  
  for (int i = 0; i<2*N/4; i++)
    {
      pa = _mm256_load_si256((__m256i *) a1); 
      pa = RED3(pa);
      _mm256_store_si256((__m256i *) a1, pa);
      a1 += 4;
    }

  

  for (int i = 0; i<2*N; i++)
    {
      //printf("%lu, ", a[i]);
      mu_assert(a[i]==a_expected[i],"Error in ntt 32");
    }
  //  printf("\n");
}


MU_TEST(ntt_64) {

  size_t N = 64;

  const integer_packed  p = BROADCAST(PRIME);
  
  const integer_packed p_64 = BROADCAST(PRIME_64);

  integer_packed _red3_t1, _red3_t2, _red3_t3, _red3_tk, retval;
 

  integer_packed pa;
   
  uint64_t *a = (uint64_t *) aligned_alloc(64, 2*N*sizeof(uint64_t));

  for (int i = 0; i<2*N; i++)
    {
      a[i] = ((uint64_t) i) * ((uint64_t) i);
    }

  /* for (int i = 0; i<2*N; i++) */
  /*   { */
  /*     printf("%lu ", a[i]); */
  /*   } */
  /* printf("\n"); */

  _Alignas(64) uint64_t a_expected[] = {85344, 62886140, 902586160, 1769346025, 605536, 1968870764, 1480421974, 1191502019, 2147481631, 759206148, 60461929, 616700930, 2147477535, 2074422418, 1147111716, 1677526598, 4160, 1408728217, 1314101250, 1536409771, 2147475519, 382548831, 1089002736, 828354172, 2147475519, 2064164197, 1548074340, 976459623, 2147479615, 2016577184, 108638100, 1201549945, 406851712, 919047998, 1889303145, 244722777, 2017452159, 966626819, 225219040, 174624300, 1740640383, 2016948929, 58867237, 281622216, 130015360, 895636476, 1818076217, 446878335, 398450816, 584553050, 1603254362, 550203204, 2009067647, 1416812235, 1437990051, 938469679, 1749016703, 774385573, 1213323029, 466924955, 138408064, 1015876404, 1283413314, 2131000868, 1532920152, 1930761614, 574574916, 274212022, 1409605185, 1742460935, 762867403, 397518797, 796704802, 1953934346, 1372920540, 1641449015, 2075691397, 134931948, 1044431099, 362454228, 1378196102, 48166360, 1738938739, 1976206587, 415636090, 1204670151, 1078479313, 987822456, 1567378358, 748255744, 709186367, 2081845875, 408270650, 1408706978, 1700158031, 896914048, 573263620, 96384392, 179759798, 1039960684, 1679675454, 1846208183, 1193241045, 396978633, 50536555, 1743486800, 1529721147, 822379771, 641411340, 960462093, 1292772850, 49250624, 736434018, 100563022, 1642697707, 2010596700, 1509554973, 337503127, 278230575, 1419713723, 1149235575, 626245989, 1647384399, 85403868, 1792209434, 1760240197, 2045036063, 1126467489};
  
  cs64(a);

  
  /* for (int i = 0; i<2*N; i++) */
  /*   { */
  /*     printf("%lu ", a[i]); */
  /*   } */
  /* printf("\n"); */
  
  uint64_t * a1 = a;
  
  for (int i = 0; i<2*N/4; i++)
    {
      pa = _mm256_load_si256((__m256i *) a1); 
      pa = RED3(pa);
      _mm256_store_si256((__m256i *) a1, pa);
      a1 += 4;
    }

  

  for (int i = 0; i<2*N; i++)
    {
      mu_assert(a[i]==a_expected[i],"Error in ntt 64");
      //      printf("%lu, ", a[i]);
    }
  //  printf("\n");

 
}


MU_TEST(ntt_128) {

  size_t N = 128;
  
  const integer_packed  p = BROADCAST(PRIME);
  
  const integer_packed p_64 = BROADCAST(PRIME_64);

  integer_packed _red3_t1, _red3_t2, _red3_t3, _red3_tk, retval;


  integer_packed pa;
   
  uint64_t *a = (uint64_t *) aligned_alloc(64, 2*N*sizeof(uint64_t));

  for (int i = 0; i<2*N; i++)
    {
      a[i] = i;
      //      printf("%lu ", a[i]);
    }
  //  printf("\n");

  //when using  cspass_sr - only the order may differ.
  _Alignas(64) uint64_t a1_expected[] = {8128, 253544855, 1505435001, 642048518, 24512, 1893938664, 642048518, 1505435001, 2147483583, 1902327400, 1149138356, 998345163, 2147483583, 245156119, 998345163, 1149138356, 0, 981095458, 303748374, 794060740, 2147483519, 1166388061, 1843735145, 1353422779, 2147483519, 1157999581, 1353422779, 1843735145, 0, 989483938, 794060740, 303748374, 4194304, 989483938, 163995505, 755522268, 2143289215, 1157999581, 1983488014, 1391961251, 2143289343, 1166388061, 1798195411, 1223445608, 4194176, 981095458, 349288108, 924037911, 4194176, 245156119, 924037911, 349288108, 2143289343, 1902327400, 1223445608, 1798195411, 2143289215, 1893938664, 1391961251, 1983488014, 4194304, 253544855, 755522268, 163995505, 386839270, 118911515, 2028572004, 1760644249, 1760644249, 2028572004, 118911515, 386839270, 476547085, 517527340, 1629956179, 1670936434, 1670936434, 1629956179, 517527340, 476547085, 208846496, 276177390, 575585087, 58053303, 1938637023, 1871306129, 1571898432, 2089430216, 2089430216, 1571898432, 1871306129, 1938637023, 58053303, 575585087, 276177390, 208846496, 1165761181, 1345283901, 1619472033, 592225761, 981722338, 802199618, 528011486, 1555257758, 1589219214, 1130190756, 1226587830, 995895719, 558264305, 1017292763, 920895689, 1151587800, 1151587800, 920895689, 1017292763, 558264305, 995895719, 1226587830, 1130190756, 1589219214, 1555257758, 528011486, 802199618, 981722338, 592225761, 1619472033, 1345283901, 1165761181, 1450969439, 1787271845, 120197299, 692390971, 696514080, 360211674, 2027286220, 1455092548, 1470192748, 777904794, 63841416, 270640898, 677290771, 1369578725, 2083642103, 1876842621, 1902981633, 882110129, 867215535, 125201308, 244501886, 1265373390, 1280267984, 2022282211, 1197596184, 1149266656, 163739246, 30490645, 949887335, 998216863, 1983744273, 2116992874, 782352893, 1202351117, 1644944672, 827000906, 1365130626, 945132402, 502538847, 1320482613, 1908214909, 2001155502, 1054893755, 1558305771, 239268610, 146328017, 1092589764, 589177748, 718944714, 1263347570, 1814151585, 99107981, 1428538805, 884135949, 333331934, 2048375538, 1541436798, 578443808, 1329645279, 935946699, 606046721, 1569039711, 817838240, 1211536820, 696514080, 998216863, 1455092548, 2027286220, 1450969439, 1149266656, 692390971, 120197299, 677290771, 1265373390, 1876842621, 2083642103, 1470192748, 882110129, 270640898, 63841416, 949887335, 1369578725, 2116992874, 1983744273, 1197596184, 777904794, 30490645, 163739246, 244501886, 360211674, 2022282211, 1280267984, 1902981633, 1787271845, 125201308, 867215535, 1320482613, 502538847, 945132402, 1365130626, 827000906, 1644944672, 1202351117, 782352893, 589177748, 1092589764, 146328017, 239268610, 1558305771, 1054893755, 2001155502, 1908214909, 1211536820, 817838240, 1569039711, 606046721, 935946699, 1329645279, 578443808, 1541436798, 2048375538, 333331934, 884135949, 1428538805, 99107981, 1814151585, 1263347570, 718944714};


  //when using cspass_rad8 - only the order may differ from previous values
  _Alignas(64) uint64_t a2_expected[] = {8128, 1505435001, 386839270, 1760644249, 24512, 642048518, 1760644249, 386839270, 2147483583, 1149138356, 476547085, 1670936434, 2147483583, 998345163, 1670936434, 476547085, 0, 303748374, 208846496, 58053303, 2147483519, 1843735145, 1938637023, 2089430216, 2147483519, 1353422779, 2089430216, 1938637023, 0, 794060740, 58053303, 208846496, 4194304, 163995505, 1165761181, 592225761, 2143289215, 1983488014, 981722338, 1555257758, 2143289343, 1798195411, 1589219214, 995895719, 4194176, 349288108, 558264305, 1151587800, 2143289215, 1391961251, 1555257758, 981722338, 4194304, 755522268, 592225761, 1165761181, 4194176, 924037911, 1151587800, 558264305, 2143289343, 1223445608, 995895719, 1589219214, 253544855, 755522268, 1345283901, 2028572004, 1893938664, 1391961251, 802199618, 118911515, 1902327400, 1223445608, 1130190756, 1629956179, 245156119, 924037911, 1017292763, 517527340, 1157999581, 1983488014, 528011486, 1871306129, 989483938, 163995505, 1619472033, 276177390, 981095458, 349288108, 920895689, 575585087, 1166388061, 1798195411, 1226587830, 1571898432, 989483938, 794060740, 276177390, 1619472033, 1157999581, 1353422779, 1871306129, 528011486, 1166388061, 1843735145, 1571898432, 1226587830, 981095458, 303748374, 575585087, 920895689, 1893938664, 642048518, 118911515, 802199618, 253544855, 1505435001, 2028572004, 1345283901, 245156119, 998345163, 517527340, 1017292763, 1902327400, 1149138356, 1629956179, 1130190756, 1450969439, 827000906, 1320482613, 696514080, 696514080, 1320482613, 827000906, 1450969439, 1470192748, 1558305771, 589177748, 677290771, 677290771, 589177748, 1558305771, 1470192748, 1902981633, 99107981, 1211536820, 949887335, 244501886, 2048375538, 935946699, 1197596184, 1197596184, 935946699, 2048375538, 244501886, 949887335, 1211536820, 99107981, 1902981633, 1787271845, 1644944672, 817838240, 998216863, 360211674, 502538847, 1329645279, 1149266656, 777904794, 1054893755, 333331934, 1265373390, 1369578725, 1092589764, 1814151585, 882110129, 1149266656, 1329645279, 502538847, 360211674, 998216863, 817838240, 1644944672, 1787271845, 882110129, 1814151585, 1092589764, 1369578725, 1265373390, 333331934, 1054893755, 777904794, 120197299, 782352893, 945132402, 1455092548, 2027286220, 1365130626, 1202351117, 692390971, 63841416, 1908214909, 146328017, 1876842621, 2083642103, 239268610, 2001155502, 270640898, 163739246, 1541436798, 884135949, 2022282211, 1983744273, 606046721, 1263347570, 125201308, 867215535, 718944714, 1569039711, 2116992874, 1280267984, 1428538805, 578443808, 30490645, 125201308, 1263347570, 606046721, 1983744273, 2022282211, 884135949, 1541436798, 163739246, 30490645, 578443808, 1428538805, 1280267984, 2116992874, 1569039711, 718944714, 867215535, 692390971, 1202351117, 1365130626, 2027286220, 1455092548, 945132402, 782352893, 120197299, 270640898, 2001155502, 239268610, 2083642103, 1876842621, 146328017, 1908214909, 63841416};

   
  cs128(a);

  
  /* for (int i = 0; i<2*N; i++) */
  /*   { */
  /*     printf("%lu ", a[i]); */
  /*   } */
  /* printf("\n"); */
  
  uint64_t * a1 = a;
  
  for (int i = 0; i<2*N/4; i++)
    {
      pa = _mm256_load_si256((__m256i *) a1); 
      pa = RED3(pa);
      _mm256_store_si256((__m256i *) a1, pa);
      a1 += 4;
    }


  for (int i = 0; i<2*N; i++)
    {
      //printf("%lu, ", a[i]);
      mu_assert((a[i]==a1_expected[i] || a[i]==a2_expected[i]) ,"Error in ntt 128");
    }
  //  printf("\n");


}



MU_TEST(inv_ntt_64) {

  size_t N = 64;

  const integer_packed  p = BROADCAST(PRIME);
  
  const integer_packed p_64 = BROADCAST(PRIME_64);

  integer_packed _mul1_t1, _mul1_t2, _mul1_t3;
 
  integer_packed _red3_t1, _red3_t2, _red3_t3, _red3_tk, retval;
 

  integer_packed pa;
   
  uint64_t *a = (uint64_t *) aligned_alloc(64, 2*N*sizeof(uint64_t));

  for (int i = 0; i<2*N; i++)
    {
      a[i] = ((uint64_t) i) * ((uint64_t) i);
    }

  /* for (int i = 0; i<2*N; i++) */
  /*   { */
  /*     printf("%lu ", a[i]); */
  /*   } */
  /* printf("\n"); */

  _Alignas(64) uint64_t a_expected[] = {85344, 62886140, 902586160, 1769346025, 605536, 1968870764, 1480421974, 1191502019, 2147481631, 759206148, 60461929, 616700930, 2147477535, 2074422418, 1147111716, 1677526598, 4160, 1408728217, 1314101250, 1536409771, 2147475519, 382548831, 1089002736, 828354172, 2147475519, 2064164197, 1548074340, 976459623, 2147479615, 2016577184, 108638100, 1201549945, 406851712, 919047998, 1889303145, 244722777, 2017452159, 966626819, 225219040, 174624300, 1740640383, 2016948929, 58867237, 281622216, 130015360, 895636476, 1818076217, 446878335, 398450816, 584553050, 1603254362, 550203204, 2009067647, 1416812235, 1437990051, 938469679, 1749016703, 774385573, 1213323029, 466924955, 138408064, 1015876404, 1283413314, 2131000868, 1532920152, 1930761614, 574574916, 274212022, 1409605185, 1742460935, 762867403, 397518797, 796704802, 1953934346, 1372920540, 1641449015, 2075691397, 134931948, 1044431099, 362454228, 1378196102, 48166360, 1738938739, 1976206587, 415636090, 1204670151, 1078479313, 987822456, 1567378358, 748255744, 709186367, 2081845875, 408270650, 1408706978, 1700158031, 896914048, 573263620, 96384392, 179759798, 1039960684, 1679675454, 1846208183, 1193241045, 396978633, 50536555, 1743486800, 1529721147, 822379771, 641411340, 960462093, 1292772850, 49250624, 736434018, 100563022, 1642697707, 2010596700, 1509554973, 337503127, 278230575, 1419713723, 1149235575, 626245989, 1647384399, 85403868, 1792209434, 1760240197, 2045036063, 1126467489};
  
  cs64(a);

  
  /* for (int i = 0; i<2*N; i++) */
  /*   { */
  /*     printf("%lu ", a[i]); */
  /*   } */
  /* printf("\n"); */
  
  uint64_t * a1 = a;
  
  for (int i = 0; i<2*N/4; i++)
    {
      pa = _mm256_load_si256((__m256i *) a1); 
      pa = RED3(pa);
      _mm256_store_si256((__m256i *) a1, pa);
      a1 += 4;
    }

  

  for (int i = 0; i<2*N; i++)
    {
      mu_assert(a[i]==a_expected[i],"Error in ntt 64");
      //      printf("%lu, ", a[i]);
    }
  //  printf("\n");

 

  us64(a);

  integer_packed inv_N = BROADCAST(33554432UL);

  
  a1 = a;
  
  for (int i = 0; i<2*N/4; i++)
    {
      pa = _mm256_load_si256((__m256i *) a1);
      pa = MUL1(pa, inv_N);
      pa = RED3(pa);
      _mm256_store_si256((__m256i *) a1, pa);
      a1 += 4;
    }


  for (int i = 0; i<2*N; i++)
    {
      mu_assert(a[i]==((uint64_t) i) * ((uint64_t) i),"Error in ntt 64");
      //printf("%lu, ", a[i]);
    }
  //  printf("\n");

  
 
}


MU_TEST(inv_ntt_128) {

  size_t N = 128;
  
  const integer_packed  p = BROADCAST(PRIME);
  
  const integer_packed p_64 = BROADCAST(PRIME_64);

  integer_packed _mul1_t1, _mul1_t2, _mul1_t3;

  integer_packed _red3_t1, _red3_t2, _red3_t3, _red3_tk, retval;


  integer_packed pa;
   
  uint64_t *a = (uint64_t *) aligned_alloc(64, 2*N*sizeof(uint64_t));

  for (int i = 0; i<2*N; i++)
    {
      a[i] = i;
      //      printf("%lu ", a[i]);
    }
  //  printf("\n");

 
   
  cs128(a);

  
  /* for (int i = 0; i<2*N; i++) */
  /*   { */
  /*     printf("%lu ", a[i]); */
  /*   } */
  /* printf("\n"); */
  
  uint64_t * a1 = a;
  
  for (int i = 0; i<2*N/4; i++)
    {
      pa = _mm256_load_si256((__m256i *) a1); 
      pa = RED3(pa);
      _mm256_store_si256((__m256i *) a1, pa);
      a1 += 4;
    }


  us128(a);

  integer_packed inv_N = BROADCAST(16777216UL);

  
  a1 = a;
  
  for (int i = 0; i<2*N/4; i++)
    {
      pa = _mm256_load_si256((__m256i *) a1);
      pa = MUL1(pa, inv_N);
      pa = RED3(pa);
      _mm256_store_si256((__m256i *) a1, pa);
      a1 += 4;
    }


  for (int i = 0; i<2*N; i++)
    {
      mu_assert(a[i]==(uint64_t) i,"Error in ntt 64");
      // printf("%lu, ", a[i]);
    }
  //  printf("\n");

 

}



MU_TEST(real_ntt_128) {

  size_t N = 128;
  
  const integer_packed  p = BROADCAST(PRIME);
  
  const integer_packed p_64 = BROADCAST(PRIME_64);

  integer_packed _red3_t1, _red3_t2, _red3_t3, _red3_tk, retval;


  integer_packed pa;
   
  uint64_t *a = (uint64_t *) aligned_alloc(64, N*sizeof(uint64_t));

  for (int i = 0; i<N; i++)
    {
      a[i] = i * i;
      //   printf("%lu ", a[i]);
    }
  //  printf("\n");

  r128(a);


  uint64_t * a1 = a;
  
  for (int i = 0; i<N/4; i++)
    {
      pa = _mm256_load_si256((__m256i *) a1); 
      pa = RED3(pa);
      _mm256_store_si256((__m256i *) a1, pa);
      a1 += 4;
    }

  _Alignas(64) uint64_t a_expected[] ={690880, 837267617, 524444442, 1458041710, 2147475519, 1905988720, 577839910, 2077389266, 2147475583, 1343755102, 677158763, 1884872482, 8192, 1315220366, 1086653883, 165260215, 2139087103, 1278189403, 143073130, 609066345, 536862720, 1121308357, 1922389229, 388270571, 8380672, 835724452, 936996368, 208740357, 536879104, 2099900730, 708051503, 1664080009, 1807128078, 658663392, 1743764695, 357848787, 2024172776, 1106415930, 1749827887, 1959187064, 290666074, 2061129159, 1238696764, 1179367143, 1278990691, 590878881, 1364463036, 328485345, 1312554234, 1558809693, 1630280073, 757348631, 1184927731, 773125051, 236944201, 1156507887, 1396097202, 41708516, 806001691, 339706740, 988380035, 642977955, 1133998304, 660455330, 1132670743, 484408868, 762111993, 783859921, 1106980905, 1009945877, 1794606712, 1568260694, 1800890659, 470857158, 1872557959, 1530474174, 793881000, 1359909585, 418225148, 1865178863, 1215220409, 2065002415, 572246374, 1452450905, 1231478542, 906528587, 665552972, 1154093560, 2094943918, 1371919130, 1229835705, 487836957, 1326502840, 1070231483, 516204790, 392156542, 1403019875, 670412145, 1837675592, 562314890, 629425457, 543395791, 1429612231, 1356930734, 28391657, 688431199, 1386347556, 2117792980, 1894651356, 283177196, 597358455, 1585864433, 336457218, 172892125, 1865589923, 375539688, 1689585079, 1882422718, 1026371618, 1882906459, 86060444, 709035136, 360579459, 284296458, 1948411525, 99528802, 647223627, 1830618410};

  for (int i = 0; i<N; i++)
    {
      //     printf("%lu, ", a[i]);
      mu_assert(a[i]==a_expected[i],"Error in ntt real 128");
    }
  // printf("\n");

  
  /* for (int i = 0; i<N/8; i++) */
  /*   { */
  /*     printf("%lu %lu %lu %lu ", a[8*i+0],  a[8*i+1],  a[8*i+2],  a[8*i+3]); */
  /*   } */
  /* printf("\n"); */

  /* for (int i = 0; i<N/8; i++) */
  /*   { */
  /*     printf("%lu %lu %lu %lu ", a[8*i+4],  a[8*i+5],  a[8*i+6],  a[8*i+7]); */
  /*   } */
  /* printf("\n"); */


}


MU_TEST(real_ntt_32) {

  size_t N = 32;
  
  const integer_packed  p = BROADCAST(PRIME);
  
  const integer_packed p_64 = BROADCAST(PRIME_64);

  integer_packed _red3_t1, _red3_t2, _red3_t3, _red3_tk, retval;


  integer_packed pa;
   
  uint64_t *a = (uint64_t *) aligned_alloc(64, N*sizeof(uint64_t));

  for (int i = 0; i<N; i++)
    {
      a[i] = i*i;
      //    printf("%lu ", a[i]);
    }
  // printf("\n");

  r32(a);


  uint64_t * a1 = a;
  
  for (int i = 0; i<N/4; i++)
    {
      pa = _mm256_load_si256((__m256i *) a1); 
      pa = RED3(pa);
      _mm256_store_si256((__m256i *) a1, pa);
      a1 += 4;
    }

  _Alignas(64) uint64_t a_expected[] = {10416, 746189352, 1204854470, 1438253787, 2147483151, 119124295, 841421362, 398272285, 2147483167, 1409682135, 1779903962, 1544961480, 512, 1961249464, 1544310875, 949852859, 2145386047, 1930161622, 1109511642, 689139034, 33553920, 741170412, 1864979790, 1500661918, 2096704, 208932649, 234250628, 589057537, 33554944, 1473421075, 2057519138, 1311964552};
  for (int i = 0; i<N; i++)
    {
      // printf("%lu, ", a[i]);
            mu_assert(a[i]==a_expected[i],"Error in ntt real 32");
    }
  //  printf("\n");


  
  /* for (int i = 0; i<N/8; i++) */
  /*   { */
  /*     printf("%lu %lu %lu %lu ", a[8*i+0],  a[8*i+1],  a[8*i+2],  a[8*i+3]); */
  /*   } */
  /* printf("\n"); */

  /* for (int i = 0; i<N/8; i++) */
  /*   { */
  /*     printf("%lu %lu %lu %lu ", a[8*i+4],  a[8*i+5],  a[8*i+6],  a[8*i+7]); */
  /*   } */
  /* printf("\n"); */

}

MU_TEST(real_ntt_64) {

  size_t N = 64;
  
  const integer_packed  p = BROADCAST(PRIME);
  
  const integer_packed p_64 = BROADCAST(PRIME_64);

  integer_packed _red3_t1, _red3_t2, _red3_t3, _red3_tk, retval;


  integer_packed pa;
   
  uint64_t *a = (uint64_t *) aligned_alloc(64, N*sizeof(uint64_t));

  for (int i = 0; i<N; i++)
    {
      a[i] = i;
      //      printf("%lu ", a[i]);
    }
  //  printf("\n");

  r64(a);


  uint64_t * a1 = a;
  
  for (int i = 0; i<N/4; i++)
    {
      pa = _mm256_load_si256((__m256i *) a1); 
      pa = RED3(pa);
      _mm256_store_si256((__m256i *) a1, pa);
      a1 += 4;
    }

  _Alignas(64) uint64_t a_expected[] = {2016, 2147483615, 2147483615, 2147483615, 2147483615, 946969364, 321024291, 1769722481, 2147483615, 2147483615, 2147483615, 2147483615, 32, 1196319915, 1572914437, 1535760811, 2147483615, 2147483615, 2147483615, 2147483615, 2097120, 1656935886, 1995609428, 1972839561, 2147483615, 2147483615, 2147483615, 2147483615, 2097184, 494742001, 397030402, 1155739608, 2147483615, 2147483615, 2147483615, 2147483615, 1954063980, 490861201, 401099841, 1014286034, 2147483615, 2147483615, 2147483615, 2147483615, 835468249, 1352874008, 1582388237, 1888719945, 2147483615, 2147483615, 2147483615, 2147483615, 2043060367, 1571689715, 613293947, 2009394920, 2147483615, 2147483615, 2147483615, 2147483615, 1102768507, 1369854736, 1883477872, 1361534399};
    
    for (int i = 0; i<N; i++)
      {
	//printf("%lu, ", a[i]);
      mu_assert(a[i]==a_expected[i],"Error in ntt real 32");
    }
    //    printf("\n");

    
  /* for (int i = 0; i<N/8; i++) */
  /*   { */
  /*     printf("%lu %lu %lu %lu ", a[8*i+0],  a[8*i+1],  a[8*i+2],  a[8*i+3]); */
  /*   } */
  /* printf("\n"); */

  /* for (int i = 0; i<N/8; i++) */
  /*   { */
  /*     printf("%lu %lu %lu %lu ", a[8*i+4],  a[8*i+5],  a[8*i+6],  a[8*i+7]); */
  /*   } */
  /* printf("\n"); */

}


MU_TEST(real_inv_ntt_128) {

  size_t N = 128;
  
  const integer_packed  p = BROADCAST(PRIME);
  
  const integer_packed p_64 = BROADCAST(PRIME_64);

  integer_packed _mul1_t1, _mul1_t2, _mul1_t3;
    
  integer_packed _red3_t1, _red3_t2, _red3_t3, _red3_tk, retval;


  integer_packed pa;
   
  uint64_t *a = (uint64_t *) aligned_alloc(64, N*sizeof(uint64_t));

  for (int i = 0; i<N; i++)
    {
      a[i] = i;
      //   printf("%lu ", a[i]);
    }
  //  printf("\n");

  r128(a);


  uint64_t * a1 = a;
  
  for (int i = 0; i<N/4; i++)
    {
      pa = _mm256_load_si256((__m256i *) a1); 
      pa = RED3(pa);
      _mm256_store_si256((__m256i *) a1, pa);
      a1 += 4;
    }

  /* for (int i = 0; i<N; i++) */
  /*   { */
  /*     //    a[i] = i * i; */
  /*     printf("%lu ", a[i]); */
  /*   } */
  /* printf("\n"); */

  v128(a);

  integer_packed inv_N = BROADCAST(33554432UL);

  for (int i = 0; i<N; i++)
    {
      // printf("%lu ", a[i]);
    }
  // printf("\n");
  
  
  a1 = a;
  
  for (int i = 0; i<N/4; i++)
    {
      pa = _mm256_load_si256((__m256i *) a1);
      pa = MUL1(pa, inv_N);
      pa = RED3(pa);
      _mm256_store_si256((__m256i *) a1, pa);
      a1 += 4;
    }


  for (int i = 0; i<N; i++)
    {
      mu_assert(a[i]==(uint64_t) i,"Error in ntt 64");
      //printf("%lu, ", a[i]);
    }
  //  printf("\n");


  
}


MU_TEST(real_inv_ntt_32) {

  size_t N = 32;
  
  const integer_packed  p = BROADCAST(PRIME);
  
  const integer_packed p_64 = BROADCAST(PRIME_64);

  integer_packed _mul1_t1, _mul1_t2, _mul1_t3;
    
  integer_packed _red3_t1, _red3_t2, _red3_t3, _red3_tk, retval;


  integer_packed pa;
   
  uint64_t *a = (uint64_t *) aligned_alloc(64, N*sizeof(uint64_t));

  for (int i = 0; i<N; i++)
    {
      a[i] = i;
      //    printf("%lu ", a[i]);
    }
  //   printf("\n");

  r32(a);


  uint64_t * a1 = a;
  
  for (int i = 0; i<N/4; i++)
    {
      pa = _mm256_load_si256((__m256i *) a1);
      pa = RED3(pa);
      _mm256_store_si256((__m256i *) a1, pa);
      a1 += 4;
    }

  /* for (int i = 0; i<N; i++) */
  /*   { */
  /*     //    a[i] = i * i; */
  /*     printf("%lu ", a[i]); */
  /*   } */
  /* printf("\n"); */

  v32(a);

  integer_packed inv_N = BROADCAST(134217728UL);

  for (int i = 0; i<N; i++)
    {
      
      //  printf("%lu ", a[i]);
    }
  //  printf("\n");
  
  
  a1 = a;
  
  for (int i = 0; i<N/4; i++)
    {
      pa = _mm256_load_si256((__m256i *) a1);
      pa = MUL1(pa, inv_N);
      pa = RED3(pa);
      _mm256_store_si256((__m256i *) a1, pa);
      a1 += 4;
    }


  for (int i = 0; i<N; i++)
    {
      //  printf("%lu, ", a[i]);
      mu_assert(a[i]==(uint64_t) i,"Error in inverse real ntt 32");
      
    }
  // printf("\n");


  
}



MU_TEST(test_real_untwisted) {

  int i;

  uint64_t dAr[4] = {0,1,4,9};
  uint64_t dBr[4] = {10,15,26,37};

  uint64_t a2[4];
  uint64_t a3[4];

  uint64_t res_Ar[4];
  uint64_t res_Br[4];
  
  integer_packed Ar, Br, c, s;

  c = LOAD(d32[2*V]);
  s = LOAD(d32[3*V]);
  
  Ar = LOAD(dAr[0]);
  Br = LOAD(dBr[0]);

  real_twisted(Ar,Br,c,s,a2,a3);

  /* for (i=0; i<4; i++) */
  /*   { */
  /*     printf("%lu %lu  ", a2[i], a3[i]); */
  /*   } */

  /* printf("\n"); */

  real_untwisted(&Ar,&Br,c,s,a2,a3);

  STORE(res_Ar[0], Ar);
  STORE(res_Br[0], Br);

  for (i=0; i<4; i++)
    {
      //      printf("%lu %lu  ", res_Ar[i]-PRIME, dAr[i]);
      mu_assert(dAr[i] == res_Ar[i]-PRIME, "Error in untwisted");
      mu_assert(dBr[i] == res_Br[i]-PRIME,  "Error in untwisted");
      
    }
  //  printf("\n");
  
  
}


MU_TEST(test_vpass_sr) {

  const integer_packed  p = BROADCAST(PRIME);
  
  const integer_packed p_64 = BROADCAST(PRIME_64);
  
  // integer_packed _mul1_t1, _mul1_t2, _mul1_t3;
 
  integer_packed _red3_t1, _red3_t2, _red3_t3, _red3_tk, retval;
 

  integer_packed pa;
   
  int i;

  int N = 32;

  uint64_t *a = aligned_alloc(32, N*sizeof(uint64_t));
  uint64_t *b = aligned_alloc(32, N*sizeof(uint64_t));

  for (i=0; i<N; i++)
    {
      a[i] = i*i;
      b[i] = a[i];
    }

  rpass_sr(a, d32, 4);

 
  //No need to multiply first part of array by 2, however the final
  //result is only multiplied by 4 instead of 8 for complex case
 
  vpass_sr(a, d32, 4);

  uint64_t * a1 = a;
  
  for (int i = 0; i<N/4; i++)
    {
      pa = _mm256_load_si256((__m256i *) a1); 
      pa = RED3(pa);
      _mm256_store_si256((__m256i *) a1, pa);
      a1 += 4;
    }

  
  for (i=0; i<N; i++)
    {
      //    printf("%lu ", a[i]);
      mu_assert(a[i] == 2*b[i] ,"Error in vpass_sr");
    }
  //  printf("\n");

}



MU_TEST(test_vpass_srl3) {

  const integer_packed  p = BROADCAST(PRIME);
  
  const integer_packed p_64 = BROADCAST(PRIME_64);
  
  //integer_packed _mul1_t1, _mul1_t2, _mul1_t3;
 
  integer_packed _red3_t1, _red3_t2, _red3_t3, _red3_tk, retval;
 

  integer_packed pa;
   
  int i;

  int N = 32;

  uint64_t *a = aligned_alloc(32, N*sizeof(uint64_t));
  uint64_t *b = aligned_alloc(32, N*sizeof(uint64_t));

  for (i=0; i<N; i++)
    {
      a[i] = i;
      b[i] = a[i];
    }

  rpass_srl3(a, d16, d32, 4);
 
 
  //No need to multiply first part of array by 2, however the final
  //result is only multiplied by 4 instead of 8 for complex case

  vpass_srl3(a, d16, d32, 4);

  uint64_t * a1 = a;
  
  for (int i = 0; i<N/4; i++)
    {
      pa = _mm256_load_si256((__m256i *) a1); 
      pa = RED3(pa);
      _mm256_store_si256((__m256i *) a1, pa);
      a1 += 4;
    }

  
  for (i=0; i<N; i++)
    {
      // printf("%lu ", a[i]);
      mu_assert(a[i] == 4*b[i] ,"Error in vpass_sr");
    }
  // printf("\n");

}


void finalize_r8_from_c4(integer * x0,  integer * y0, integer *  x1, integer * y1, integer *  x2, integer *   y2, integer *   x3, integer *  y3);

void unfinalize_c4_from_r8(integer * a0,  integer * b0, integer *  a1, integer * b1, integer *  a2, integer *   b2, integer *   a3, integer *  b3);


MU_TEST(test_unfinalize_c4_from_r8)
{
  integer x0,x1,x2,x3,y0,y1,y2,y3;

  x0 = 1*1;
  y0 = 2*2;

  x1 = 3*3;
  y1 = 4*4;

  x2 = 5*5;
  y2 = 6*6;

  x3 = 7*7;
  y3 = 8*8;

  finalize_r8_from_c4(&x0, &y0, &x1, &y1, &x2, &y2, &x3, &y3);

  unfinalize_c4_from_r8(&x0, &y0, &x1, &y1, &x2, &y2, &x3, &y3);


  //  printf("%lu %lu %lu %lu %lu %lu %lu %lu\n", x0-PRIME,y0-PRIME,x1-PRIME,y1-PRIME,x2,y2-PRIME,x3-PRIME,y3-PRIME);
  /* mu_assert_double_close_rel(x0, 1*1, 4e-16); */
  /* mu_assert_double_close_rel(y0, 2*2, 4e-16); */
  /* mu_assert_double_close_rel(x1, 3*3, 4e-16); */
  /* mu_assert_double_close_rel(y1, 4*4, 5e-16); */
  /* mu_assert_double_close_rel(x2, 5*5, 4e-16); */
  /* mu_assert_double_close_rel(y2, 6*6, 4e-16); */
  /* mu_assert_double_close_rel(x3, 7*7, 4e-16); */
  /* mu_assert_double_close_rel(y3, 8*8, 4e-16); */
  
  
}

MU_TEST(test_scalec)
{
  const integer_packed  p = BROADCAST(PRIME);
  const integer_packed p_64 = BROADCAST(PRIME_64);
  
  const size_t n = 32;
  integer a[2*n];

  uint64_t * a1 = a;

  integer_packed pa;
  integer_packed _red3_t1, _red3_t2, _red3_t3, _red3_tk, retval;
  
  int i;

  for (i=0; i<2*n; i++)
    {
      a[i] = i*32;
    }

  /* for (i=0; i<2*n; i++) */
  /*   { */
  /*     printf("%lu ", a[i]); */
  /*   } */
  /* printf("\n"); */
  
       
  scalec(a, 32);

 
  /* for (i=0; i<2*n; i++) */
  /*   { */
  /*     printf("%lu ", a[i]); */
  /*   } */
  /* printf("\n"); */
  
  
  for (int i = 0; i<2*n/4; i++)
    {
      pa = _mm256_load_si256((__m256i *) a1); 
      pa = RED3(pa);
      _mm256_store_si256((__m256i *) a1, pa);
      a1 += V;
    }


  for (i=0; i<2*n; i++)
    {
      mu_assert(a[i] == i, "Error in scalec");
      // printf("%lu ", a[i]);
    }
  //  printf("\n");
  
}


MU_TEST(test_scaler)
{
  const integer_packed  p = BROADCAST(PRIME);
  const integer_packed p_64 = BROADCAST(PRIME_64);
  
  const size_t n = 32;
  integer a[n];

  uint64_t * a1 = a;

  integer_packed pa;
  integer_packed _red3_t1, _red3_t2, _red3_t3, _red3_tk, retval;
  
  int i;

  for (i=0; i<n; i++)
    {
      a[i] = (i+1)*32;
    }

  /* for (i=0; i<n; i++) */
  /*   { */
  /*     printf("%lu ", a[i]); */
  /*   } */
  /* printf("\n"); */
  
       
  scaler(a, 32);

 
  /* for (i=0; i<n; i++) */
  /*   { */
  /*     printf("%lu ", a[i]); */
  /*   } */
  /* printf("\n"); */
  
  
  for (int i = 0; i<n/4; i++)
    {
      pa = _mm256_load_si256((__m256i *) a1); 
      pa = RED3(pa);
      _mm256_store_si256((__m256i *) a1, pa);
      a1 += V;
    }


  for (i=0; i<n; i++)
    {
      //    printf("%lu ", a[i]);
      if (i==0 || i == V)
	      mu_assert(a[i] == i+1, "Error in scaler");
      else
	mu_assert(a[i] == 2*(i+1), "Error in scaler");

      
    }
  //   printf("\n");
  
}


MU_TEST_SUITE(test_suite) {
  MU_RUN_TEST(test_mul_mod_avx256);
  MU_RUN_TEST(test_reduction_avx256);
  MU_RUN_TEST(test_macros_avx256);
  MU_RUN_TEST(transpose_4x4);
  MU_RUN_TEST(rad8); 
  MU_RUN_TEST(ntt_16); 
  MU_RUN_TEST(ntt_32);
  MU_RUN_TEST(ntt_64);
  MU_RUN_TEST(ntt_128);
  MU_RUN_TEST(inv_ntt_64);
  MU_RUN_TEST(inv_ntt_128);
  MU_RUN_TEST(real_ntt_128);
  MU_RUN_TEST(real_ntt_32);
  MU_RUN_TEST(real_ntt_64);

  MU_RUN_TEST(real_inv_ntt_128);
  MU_RUN_TEST(real_inv_ntt_32);
  MU_RUN_TEST(test_real_untwisted);
  MU_RUN_TEST(test_vpass_sr);
  MU_RUN_TEST(test_vpass_srl3);
  MU_RUN_TEST(test_unfinalize_c4_from_r8);

  MU_RUN_TEST(test_scalec);
  MU_RUN_TEST(test_scaler);

}

int main(int argc, char *argv[]) {
        MU_RUN_SUITE(test_suite);
        MU_REPORT();
        return minunit_fail;
}
