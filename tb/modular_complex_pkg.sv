package modular_complex_pkg;

   typedef struct  {
      longint	   re;
      longint	   im;
  } complex_mod_t;
   
   
  function automatic complex_mod_t cmpl_add_mod(input complex_mod_t a, input complex_mod_t b, input longint p);

    complex_mod_t result;
    longint re, im;

    re = (a.re + b.re) % p;
    im = (a.im + b.im) % p;

    result.re = re;
    result.im = im;

    return result;

  endfunction


  function automatic complex_mod_t cmpl_sub_mod(input complex_mod_t a, input complex_mod_t b, input longint p);

    complex_mod_t result;
    longint re, im;


    if (a.re >= b.re) begin
      re = (a.re - b.re) % p;
    end else begin
      re = (p + a.re - b.re) % p;
    end

    if (a.im >= b.im) begin
      im = (a.im - b.im) % p;
    end else begin
      im = (p + a.im - b.im) % p;
    end

    result.re = re;
    result.im = im;
    return result;

  endfunction


  function automatic complex_mod_t cmpl_mult_mod(input complex_mod_t a, input complex_mod_t b, input longint p);
     
    complex_mod_t result;
    longint re0, re1, im0, im1, re, im;


    re0 = (a.re * b.re) % p;
    re1 = (a.im * b.im) % p;

    if (re0 >= re1) begin
      re = (re0 - re1) % p;
    end else begin
      re = (p + re0 - re1) % p;
    end
     
     im0 = (a.re * b.im) % p;
     im1 = (a.im * b.re) % p;

     im = (im0 + im1) % p;
     
     result.re = re;
     result.im = im;

    return result;

  endfunction


  function automatic complex_mod_t cmpl_power_mod(input complex_mod_t base, input integer exponent,
                                       input longint p);
    complex_mod_t result;
    complex_mod_t b;
    longint e;

    // Initialize
    result.re = 1;
    result.im = 0;
     
    b.re = base.re % p;
    b.im = base.im % p;
     
    e = exponent;

    // Modular exponentiation by squaring
    while (e > 0) begin
      if (e[0] & 1) begin
        result = cmpl_mult_mod(result, b, p);
      end
      b = cmpl_mult_mod(b, b, p);
      e = e >> 1;
    end

    return result;

  endfunction

endpackage
