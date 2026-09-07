function res = mult_mod_complex(x,y,p)

  %%  mult_mod_complex(x,y,p) - returns z = x .* y mod p
  %%
  %%  x, y, z are "complex" arrays of remainders modulo p.
  %%
  %%  first column is real part, second column is imaginary part.
  
  a0 = x(:,1);
  b0 = x(:,2);

  a1 = y(:,1);
  b1 = y(:,2);

  res = [mod(mod(a0.*a1,p) - mod(b0.*b1,p),p), mod(mod(a0.*b1,p) + mod(a1.*b0,p),p)];

end



%!test
%! p = int64(2)^13-1;
%! x = int64(floor(rand(16,2)*p));
%! y = int64(floor(rand(16,2)*p));
%! z = mult_mod_complex(x,y,p)
%! 
%! xd = double(x(:,1))+j*double(x(:,2));
%! yd = double(y(:,1))+j*double(y(:,2));
%! zd = xd.*yd;
%! z_expected = [mod(real(zd),p), mod(imag(zd),p)]
%! assert( norm(double(z-z_expected))==0)
