#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "mrsn_internal.h"

int main()
{
  int N_rand = 20000000;
  
  uint32_t a;

  uint32_t b, c, c_expected;

  uint32_t p = 2147483647u;

  uint32_t error;

  int32_t d, d_expected;
  
  //Test remainders close to 0 and p-1
  for (int i = 0; i < 1000; i++)
    {
      a = i%2 == 0 ? i/2 : p - (i-1)/2;

      c = red1(a);
      c_expected = ((uint64_t) a) % p;
      error = c - c_expected;
      if (error != 0 ) printf("a = %u, a_red = %u, expected = %u, error: %u\n", a, c, c_expected, error);

      d = red2(c);
      d_expected = ((int32_t) c) <= (((int32_t) p)-1)/2 ? ((int32_t) c) : ((int32_t) c) -  ((int32_t) p);
	
      if (d!=d_expected) printf("a = %u, a_red = %u, a_red2 = %d, expected = %d\n", a, c, d, d_expected);

      c = rot15(a);
      c_expected = ((uint64_t) a * 32768lu ) % p;
      error = c - c_expected;
      if (error != 0 && c_expected != 0 && c != p)  printf("a = %u, a*2^15 = %u, expected = %u, error: %u\n", a, c, c_expected, error);

      c = rot15_alt(a);
      c_expected = ((uint64_t) a * 32768lu ) % p;
      error = c - c_expected;
      if (error != 0 && c_expected != 0 && c != p)  printf("a = %u, a*2^15 = %u, expected = %u, error: %u\n", a, c, c_expected, error);
      

      for (int j = 0; j < 100; j++)
	{
	  b = j%2 == 0 ? j/2 : p - (j-1)/2;

	  c = add(a, b);
	  c_expected = ((uint64_t) a + (uint64_t) b) % p;
	  error = c - c_expected;
	  if (error != 0 &&  c_expected != 0 && c != p)  printf("a = %u, b = %u, a+b =` %u, expected = %u, error: %u\n", a, b, c, c_expected, error);
	  
	  c = sub(a, b);
	  c_expected = a >= b ?  ((uint64_t) a - (uint64_t) b) % p : ((uint64_t) a + ((uint64_t) p  - (uint64_t) b)) % p;
	  error = c - c_expected;
	  if (error != 0   && c_expected != 0 && c != p) printf("a = %u, b = %u, a-b = %u, expected = %u, error: %u\n", a, b, c, c_expected, error);

	  c = prod(a, b);
	  c_expected = ((uint64_t) a * (uint64_t) b) % p;
	  error = c - c_expected;
	  if (error != 0   && c_expected != 0 && c != p)  printf("a = %u, b = %u, a*b = %u, expected = %u, error: %u\n", a, b, c, c_expected, error);
	}
    }

  //Test random remainders
  for (int i = 0; i < N_rand; i++)
    {
      a = ((rand() & 0x7fffu)<<16 | (rand() & 0x7fffu)<<2 ) | (rand() & 0x7fffu)>>13;

      b = ((rand() & 0x7fffu)<<16 | (rand() & 0x7fffu)<<2 ) | (rand() & 0x7fffu)>>13;

      c = red1(a);
      c_expected = ((uint64_t) a) % p;
      error = c - c_expected;
      if (error != 0 ) printf("a = %u, a_red = %u, expected = %u, error: %u\n", a, c, c_expected, error);

      d = red2(c);
      d_expected = ((int32_t) c) <= (((int32_t) p)-1)/2 ? ((int32_t) c) : ((int32_t) c) -  ((int32_t) p);
      if (d!=d_expected) printf("a = %u, a_red = %u, a_red2 = %d, expected = %d\n", a, c, d, d_expected);

      c = rot15(a);
      c_expected = ((uint64_t) a * 32768lu ) % p;
      error = c - c_expected;
      if (error != 0 && c_expected != 0 && c != p) printf("a = %u, a*2^15 = %u, expected = %u, error: %u\n", a, c, c_expected, error);

      c = rot15_alt(a);
      c_expected = ((uint64_t) a * 32768lu ) % p;
      error = c - c_expected;
      if (error != 0 && c_expected != 0 && c != p) printf("a = %u, a*2^15 = %u, expected = %u, error: %u\n", a, c, c_expected, error);

      c = rot24(a);
      c_expected = ((uint64_t) a * 16777216lu ) % p;
      error = c - c_expected;
      if (error != 0 && c_expected != 0 && c != p) printf("a = %u, a*2^24 = %u, expected = %u, error: %u\n", a, c, c_expected, error);

      
      c = add(a, b);
      c_expected = ((uint64_t) a + (uint64_t) b) % p;
      error = c - c_expected;
      if (error != 0  && c_expected != 0 && c != p )  printf("a = %u, b= %u, a+b = %u, expected = %u, error: %u\n", a, b, c, c_expected, error);

      c = sub(a, b);
      c_expected = a >= b ?  ((uint64_t) a - (uint64_t) b) % p : ((uint64_t) a + ((uint64_t) p  - (uint64_t) b)) % p; 
      error = c - c_expected;
      if (error != 0 && c_expected != 0 && c != p) printf("a = %u, b = %u, a-b = %u, expected = %u, error: %u\n", a, b, c, c_expected, error);

      c = prod(a, b);
      c_expected = ((uint64_t) a * (uint64_t) b) % p;
      error = c - c_expected;
      if (error != 0 && c_expected != 0 && c != p) printf("a = %u, b = %u, a*b = %u, expected = %u, error: %u\n", a, b, c, c_expected, error);
    }



  // Test Barrett reduction

  int32_t half_ps = 1073741823;

  int32_t ps =  2147483647;

  int32_t q = 3329;
  int32_t Rq = 2580335;

  int32_t x, x_barrett1, x_barrett2;

  x = -17624222;

  x_barrett1 = barrett1(x, q, Rq);

  printf("x = %d, x_barrett1 = %d\n", x, x_barrett1);

  x_barrett2 = barrett2(x_barrett1, q, Rq);

  printf("x = %d, x_barrett1 = %d\n", x, x_barrett2);
  
}

  
