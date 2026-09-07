function y = intt_p13(x);

  [r, c] = size(x);

  if (r == 1)
    N = c;
    N_pow = floor(log2(c));
    xx(:,1) = x';
    xx(:,2) = int64(0);
  elseif (c == 1)
    N = r;
    N_pow = floor(log2(r));
    xx(:,1) = x';
    xx(:,2) = int64(0);
  elseif (c == 2)
    N = r;
    N_pow = floor(log2(r));
    xx = x;
  else
    disp("Array x should be either pure real of size (1,N) or (N,1) or complex of size (N,2) with real part in first column and imag part in second column");
    y = NaN;
    return;
  end
  
  if (2^N_pow ~= N || N_pow > 14)
    disp("Length of array x should be a power of 2 no larger than 2^14");
    y = NaN;
    return;
  endif

  p_pow = 13;
  
  p = int64(2)^p_pow-1;

  w0 = int64([638, 3892]);

  w = power_mod_complex(w0, 2^(p_pow+1 - N_pow), p);

  sq2 = 2^((p_pow-1)/2);
  
  assert(sum(abs([  1   0] - power_mod_complex(w,2^(N_pow),  p)))==0)
  assert(sum(abs([p-1   0] - power_mod_complex(w,2^(N_pow-1),p)))==0)
  assert(sum(abs([  0   1] - power_mod_complex(w,2^(N_pow-2),p)))==0)
  assert(sum(abs([sq2 sq2] - power_mod_complex(w,2^(N_pow-3),p)))==0)

  w_inv = [w(1) mod(-w(2),p)];
  
  roots_N(1,1:2) = [1,0];


  for i=2:N
    roots_N(i,1:2) = power_mod_complex(w_inv, i-1, p);
  end

  
  M_inv=int64(zeros(N,N,2));
  
  kk = int64([1:N]-1)';
  
  for j=1:N
    po = mod((j-1)*kk,N) + 1;
    M_inv(:,j,1) =roots_N(po,1);
    M_inv(:,j,2) =roots_N(po,2);
  end

  y = mat_vect_mod_complex(M_inv,xx,p);

  N_inv = int64(2)^(p_pow-N_pow);
  
  y = mod(y * N_inv, p);
end



%!test
%!
%! re = int64([2600, 2980, 1480, 3429, 6608, 3328, 5255, 5870, 2168, 454, 2764, 1978, 826, 551, 3783, 7693]);
%!
%! im = int64([7661, 3262, 3695, 1096, 8104, 311, 1626, 4582, 2052, 1645, 6541, 5549, 7340, 4572, 7779, 1445]);
%!
%!
%! x(:,1) = re';
%! x(:,2) = im';
%! y = ntt_p13(x);
%! z = intt_p13(y);
%! 
%! assert(sum(sum(abs(z-x)))==0)

