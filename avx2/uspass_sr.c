#include "ntt.h"

#include "pre8_V4.h"

#include "i.h"


void uspass_sr(integer *a,  const integer *w,  size_t n)
{

  integer_packed x0, x1, x2, x3, y0, y1, y2, y3;
  integer_packed c, s, t02, t13, s02, s13, Ar, Br, Ai, Bi;

  size_t k, n2;
  
  
  integer *a1;
  integer *a2;
  integer *a3;

  integer *b;
  integer *b1;
  integer *b2;
  integer *b3;

  
  n2 = n << 1; //2n
  
  a1 = a  + n2;
  a2 = a1 + n2;
  a3 = a2 + n2;

  b = a + (n<<3);
  b1 = b  + n2;
  b2 = b1 + n2;
  b3 = b2 + n2;
  
  k = n >> (VPOW-1);  //2n/V

  do {
    
    c = LOAD(w[0]);
    s = LOAD(w[V]);
    
    //    UNTWISTED(Ar, Ai, Br, Bi, c, s, b[0], b1[0], b2[0], b3[0]);
    untwisted(&Ar, &Ai, &Br, &Bi, c, s, b, b1, b2, b3);

    t02 = LOAD(a[0]);
    s02 = LOAD(a2[0]);
 
    x0 = ADD(t02, Ar);
    x2 = SUB(t02, Ar);

    y0 = ADD(s02, Ai);
    y2 = SUB(s02, Ai);

    STORE(a[0], x0);
    STORE(b[0], y0);

    t13 = LOAD(a1[0]);
    s13 = LOAD(a3[0]);
    
    x1 = ADD(t13, Br);
    x3 = SUB(t13, Br);
    
    y1 = ADD(s13, Bi);
    y3 = SUB(s13, Bi);

    STORE(a1[0], x1);
    STORE(b1[0], y1);

    STORE(a2[0], x2);
    STORE(b2[0], y2);

    STORE(a3[0], x3);
    STORE(b3[0], y3);
  
    a  += V;
    a1 += V;
    a2 += V;
    a3 += V;

    b  += V;
    b1 += V;
    b2 += V;
    b3 += V;
    
    w  += V2;
    
  }  while (k -= 1);
}
