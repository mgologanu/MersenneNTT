#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include<immintrin.h>
#include <stdalign.h>

#include <time.h>

#include "minunit.h"

#include "mrsn_avx2_internal.h"


#define MRSN_PRIME 2147483647

MU_TEST(test_add) {

  const int vn = 4;
  const int n  = 8*vn;

  __m256i va, vb, vc;
  
  uint32_t * a = aligned_alloc(32, n*sizeof(uint32_t));

  uint32_t * b = aligned_alloc(32, n*sizeof(uint32_t));
  
  uint32_t * c = aligned_alloc(32, n*sizeof(uint32_t));

  uint32_t * c_expected = aligned_alloc(32, n*sizeof(uint32_t));

  int32_t p = MRSN_PRIME;

  for (int i= 0; i<8; i++)
    {
      a[i] = (p-1)-i+1;
      b[i] = (p-1)-i+1;
    }

  for (int i= 0; i<8; i++)
    {
      a[i+8] = (p-1)-i+1;
      b[i+8] = i;
    }

  for (int i= 0; i<8; i++)
    {
      a[i+2*8] = i; 
      b[i+2*8] = (p-1)-i+1;
    }

  for (int i= 0; i<8; i++)
    {
      a[i+3*8] = i; 
      b[i+3*8] = i;
    }
    
  for (int j = 0; j < vn; j++)
    {
      va = _mm256_load_si256((__m256i *) &a[0 + j*8]);
      
      vb = _mm256_load_si256((__m256i *) &b[0 + j*8]);
  
      vc = add(va, vb);
  
      _mm256_store_si256((__m256i *) &c[0 + j*8], vc);
    }
  
  for  (int i = 0; i<n; i++)
    {
      c_expected[i] = ((int64_t) a[i] + (int64_t) b[i]) % (int64_t) p;
    }
  
  for (int i= 0; i<n; i++)
    {
      //printf("%d %d \n", c[i], c_expected[i]);
      mu_assert(c[i] == c_expected[i] || (c[i]==p && c_expected[i] == 0), "Mul gresit");
    }
  //    printf("\n");
}


MU_TEST(test_add_random) {

  const int vn = 1000;

  const int n  = 8*vn;

  __m256i va, vb, vc;
  
  uint32_t * a = aligned_alloc(32, n*sizeof(uint32_t));

  uint32_t * b = aligned_alloc(32, n*sizeof(uint32_t));
  
  uint32_t * c = aligned_alloc(32, n*sizeof(uint32_t));

  uint32_t * c_expected = aligned_alloc(32, n*sizeof(uint32_t));

  int32_t p = MRSN_PRIME;

  srand((unsigned int)time(NULL));
     
  for (int i= 0; i<n; i++)
    {
      a[i] = rand();
      b[i] = rand();
    }

  for (int j = 0; j < vn; j++)
    {
      va = _mm256_load_si256((__m256i *) &a[0 + j*8]);
      
      vb = _mm256_load_si256((__m256i *) &b[0 + j*8]);
  
      vc = add(va, vb);
  
      _mm256_store_si256((__m256i *) &c[0 + j*8], vc);
    }

  for  (int i = 0; i<n; i++)
    {
      c_expected[i] = ((int64_t) a[i] + (int64_t) b[i]) % (int64_t) p;
    }
   
  for (int i= 0; i<n; i++)
    {
      //printf("%d %d \n", c[i], c_expected[i]);
      mu_assert(c[i] == c_expected[i] || (c[i]==p && c_expected[i] == 0), "Mul gresit");
    }
  //    printf("\n");
}



MU_TEST(test_sub) {

  const int vn = 4;

  const int n  = 8*vn;

  __m256i va, vb, vc;
    
  uint32_t * a = aligned_alloc(32, n*sizeof(uint32_t));

  uint32_t * b = aligned_alloc(32, n*sizeof(uint32_t));
  
  uint32_t * c = aligned_alloc(32, n*sizeof(uint32_t));

  uint32_t * c_expected = aligned_alloc(32, n*sizeof(uint32_t));

  int32_t p = MRSN_PRIME;

  for (int i= 0; i<8; i++)
    {
      a[i] = (p-1)-i+1;
      b[i] = (p-1)-i+1;
    }

  for (int i= 0; i<8; i++)
    {
      a[i+8] = (p-1)-i+1;
      b[i+8] = i;
    }

  for (int i= 0; i<8; i++)
    {
      a[i+2*8] = i; 
      b[i+2*8] = (p-1)-i+1;
    }

  for (int i= 0; i<8; i++)
    {
      a[i+3*8] = i; 
      b[i+3*8] = i;
    }
    
  for (int j = 0; j < vn; j++)
    {
      va = _mm256_load_si256((__m256i *) &a[0 + j*8]);
      
      vb = _mm256_load_si256((__m256i *) &b[0 + j*8]);
  
      vc = sub(va, vb);
  
      _mm256_store_si256((__m256i *) &c[0 + j*8], vc);
    }

  int64_t tmp;
  int sign;

  for  (int i = 0; i<n; i++)
    {
      tmp = (int64_t) a[i] - (int64_t) b[i];
      if (tmp < 0)
	{
	  sign = 1;
	  tmp = -tmp;
	}
      else
	{
	  sign = 0;
	}
      tmp = tmp % (int64_t) p;
      if (sign && tmp != 0)
	{
	  tmp = p - tmp;
	}
      c_expected[i] = tmp;
    }


  for (int i= 0; i<n; i++)
    {
      //printf("%d %d \n", c[i], c_expected[i]);
      mu_assert(c[i] == c_expected[i] || (c[i]==p && c_expected[i] == 0), "Mul gresit");
    }
  //    printf("\n");
}


MU_TEST(test_sub_random) {

  const int vn = 1000;

  const int n  = 8*vn;
   
  __m256i va, vb, vc;
  
  uint32_t * a = aligned_alloc(32, n*sizeof(uint32_t));

  uint32_t * b = aligned_alloc(32, n*sizeof(uint32_t));
  
  uint32_t * c = aligned_alloc(32, n*sizeof(uint32_t));

  uint32_t * c_expected = aligned_alloc(32, n*sizeof(uint32_t));

  int32_t p = 2147483647;

  srand((unsigned int)time(NULL));
     
  for (int i= 0; i<n; i++)
    {
      a[i] = rand();
      b[i] = rand();
    }

  for (int j = 0; j < vn; j++)
    {
      va = _mm256_load_si256((__m256i *) &a[0 + j*8]);
      
      vb = _mm256_load_si256((__m256i *) &b[0 + j*8]);
  
      vc = sub(va, vb);
  
      _mm256_store_si256((__m256i *) &c[0 + j*8], vc);
    }

  int64_t tmp;
  int sign;

  for  (int i = 0; i<n; i++)
    {
      tmp = (int64_t) a[i] - (int64_t) b[i];
      if (tmp < 0)
	{
	  sign = 1;
	  tmp = -tmp;
	}
      else
	{
	  sign = 0;
	}
      tmp = tmp % (int64_t) p;
      if (sign && tmp != 0)
	{
	  tmp = p - tmp;
	}
      c_expected[i] = tmp;
    }

  for (int i= 0; i<n; i++)
    {
      //printf("%d %d \n", c[i], c_expected[i]);
      mu_assert(c[i] == c_expected[i] || (c[i]==p && c_expected[i] == 0), "Mul gresit");
    }
  //    printf("\n");
}


MU_TEST(test_mul) {

  const int vn = 4;
  const int n  = 8*vn;

  __m256i va, vb, vc;
  
  uint32_t * a = aligned_alloc(32, n*sizeof(uint32_t));

  uint32_t * b = aligned_alloc(32, n*sizeof(uint32_t));
  
  uint32_t * c = aligned_alloc(32, n*sizeof(uint32_t));

  uint32_t * c_expected = aligned_alloc(32, n*sizeof(uint32_t));

  int32_t p = MRSN_PRIME;

  for (int i= 0; i<8; i++)
    {
      a[i] = (p-1)-i+1;
      b[i] = (p-1)-i+1;
    }

  for (int i= 0; i<8; i++)
    {
      a[i+8] = (p-1)-i+1;
      b[i+8] = i;
    }

  for (int i= 0; i<8; i++)
    {
      a[i+2*8] = i; 
      b[i+2*8] = (p-1)-i+1;
    }

  for (int i= 0; i<8; i++)
    {
      a[i+3*8] = i; 
      b[i+3*8] = i;
    }
    
  for (int j = 0; j < vn; j++)
    {
      va = _mm256_load_si256((__m256i *) &a[0 + j*8]);
      
      vb = _mm256_load_si256((__m256i *) &b[0 + j*8]);
  
      vc = mul(va, vb);
  
      _mm256_store_si256((__m256i *) &c[0 + j*8], vc);
    }
  
  for  (int i = 0; i<n; i++)
    {
      c_expected[i] = ((int64_t) a[i] * (int64_t) b[i]) % (int64_t) p;
    }
  
  for (int i= 0; i<n; i++)
    {
      // printf("%d %d \n", c[i], c_expected[i]);
      mu_assert(c[i] == c_expected[i] || (c[i]==p && c_expected[i] == 0), "Mul gresit");
    }
  //    printf("\n");
}


MU_TEST(test_mul_random) {

  const int vn = 1000;

  const int n  = 8*vn;

  __m256i va, vb, vc;
  
  uint32_t * a = aligned_alloc(32, n*sizeof(uint32_t));

  uint32_t * b = aligned_alloc(32, n*sizeof(uint32_t));
  
  uint32_t * c = aligned_alloc(32, n*sizeof(uint32_t));

  uint32_t * c_expected = aligned_alloc(32, n*sizeof(uint32_t));

  int32_t p = MRSN_PRIME;

  srand((unsigned int)time(NULL));
     
  for (int i= 0; i<n; i++)
    {
      a[i] = rand();
      b[i] = rand();
    }

  for (int j = 0; j < vn; j++)
    {
      va = _mm256_load_si256((__m256i *) &a[0 + j*8]);
      
      vb = _mm256_load_si256((__m256i *) &b[0 + j*8]);
  
      vc = mul(va, vb);
  
      _mm256_store_si256((__m256i *) &c[0 + j*8], vc);
    }

  for  (int i = 0; i<n; i++)
    {
      c_expected[i] = ((int64_t) a[i] * (int64_t) b[i]) % (int64_t) p;
    }
   
  for (int i= 0; i<n; i++)
    {
      //printf("%d %d \n", c[i], c_expected[i]);
      mu_assert(c[i] == c_expected[i] || (c[i]==p && c_expected[i] == 0), "Mul gresit");
    }
  //    printf("\n");
}



MU_TEST_SUITE(test_suite) {
  MU_RUN_TEST(test_add);
  MU_RUN_TEST(test_add_random);
  MU_RUN_TEST(test_sub);
  MU_RUN_TEST(test_sub_random);
  MU_RUN_TEST(test_mul);
  MU_RUN_TEST(test_mul_random);
}

int main(int argc, char *argv[]) {
  MU_RUN_SUITE(test_suite);
  MU_REPORT();
  return minunit_fail;
}
