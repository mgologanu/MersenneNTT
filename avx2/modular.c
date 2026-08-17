#include<stdint.h>
#include <stdio.h>


const int64_t  p = 2147483647L; 

const uint32_t p_pow = 31;

const int64_t pm1_half =  1073741823L;



int32_t mul_mod_mersenne_32( int32_t a, int32_t b)
{

  
  //const int32_t p = 2147483647; //2^31-1

 
  int32_t r;

  int64_t T;
 
  T = ((int64_t) a)*((int64_t) b);

  r  = (int32_t) (T & p) + (int32_t) (T >> p_pow);

  return r;
}


int64_t mul_mod_mersenne_64(int64_t a, int64_t b)
{
  /* returns */
  

  int64_t T, t1, t2, t3;

  T = a*b;

  t1 = T & p;

  t2 =  (T >> p_pow);

  t3 = t1 + t2;
    
  // printf("%ld %ld %ld %ld\n", T, t1, t2, t3);
  
  return t3;

}


int64_t mul_mod_mersenne_64_2(int64_t a, int64_t b)
{
 
  int64_t T = a*b;

  int64_t t1 = (T & p) +  (T >> p_pow);

  return (t1 & p) +  (t1 >> p_pow);

}

int64_t red_mersenne_64(int64_t T)
{
  /*
  int64_t t1, t2, t3;

  t1 = T & p;

  t2 =  (T >> p_pow);

  t3 = t1 + t2;
    
  // printf("%ld %ld %ld %ld\n", T, t1, t2, t3);
  
  return t3;
  */
  return (T & p) +  (T >> p_pow);
  
}

int64_t red_mersenne_64_2(int64_t T)
{
  int64_t t1 = (T & p) +  (T >> p_pow);
  return  (t1 & p) +  (t1 >> p_pow);
}


int64_t red_mersenne_64_2_final(int64_t T)
{
  int64_t t1 = (T & p) +  (T >> p_pow);

  int64_t t2 = (t1 & p) +  (t1 >> p_pow);

  int64_t t3 = t2 - p;

  if (t2 > pm1_half)
    return t3;
  else
    return t2;  
}
