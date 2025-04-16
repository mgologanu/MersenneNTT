#include "ntt.h"
#include "pre8_V4.h"
#include "i.h"

/* 
   Split radix pass for x^8n-1 - 8n real numbers 

   a: 8n complex = 16n values = in normal order (unpacked)

   w: 2n complex numbers or 4n values = [0...4n-1] in interleaved V-packed format

   Uses 1/4 of roots (no mirroring on 45 deg). 

*/
void rpass_sr(integer *a,  const integer *w,  size_t n)
{

  integer_packed x0, x1, x2, x3;
  integer_packed c, s, t02, t13, Ar, Br;

  size_t k, n2;
  
  
  integer *a1;
  integer *a2;
  integer *a3;

  n2 = n << 1; //2n
  
  a1 = a  + n2;
  a2 = a1 + n2;
  a3 = a2 + n2;

  k = n >> (VPOW-1);  //2n/V

  do {
    
    x0 = LOAD(a[0]);
    x1 = LOAD(a1[0]);
    x2 = LOAD(a2[0]);
    x3 = LOAD(a3[0]);
    
 
    c = LOAD(w[0]);
    s = LOAD(w[V]);
    
    t02 = ADD(x0, x2);
    Ar  = SUB(x0, x2);
    
    t13 = ADD(x1, x3);
    Br  = SUB(x1, x3);
    
   
    STORE(a[0], t02);
    STORE(a1[0], t13);
 

    //REAL_TWISTED(Ar, Br, c, s, a2[0], a3[0]);
  
    real_twisted(Ar, Br, c, s, a2, a3);
		 
    a  += V;
    a1 += V;
    a2 += V;
    a3 += V;

    w  += V2;
    
  }  while (k -= 1);
}



