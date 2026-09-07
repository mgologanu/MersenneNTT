function y = nega_cyc_intt_p19(x);

  
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

  
  if (2^N_pow ~= N || N_pow > 19)
    disp("Length of array x should be a power of 2 no larger than 2^19");
    y = NaN;
    return;
  endif


  N2 = 2*N;

  N2_pow = N_pow + 1;
  
  p_pow = 19;
  
  p = int64(2)^p_pow-1;

  w0 = int64([502746   91912]);

  w = power_mod_complex(w0, 2^(p_pow+1 - N2_pow), p);
 
  sq2 = 2^((p_pow-1)/2);
  
  assert(sum(abs([  1   0] - power_mod_complex(w,2^(N2_pow),  p)))==0)
  assert(sum(abs([p-1   0] - power_mod_complex(w,2^(N2_pow-1),p)))==0)
  assert(sum(abs([  0   1] - power_mod_complex(w,2^(N2_pow-2),p)))==0)
  assert(sum(abs([sq2 sq2] - power_mod_complex(w,2^(N2_pow-3),p)))==0)
  
  

  %% First method - double the original array [x -x]
  
  xx_double(1:2:N2,1:2) = zeros(N,2);
  
  xx_double(2:2:N2, 1:2) = xx(1:N,1:2);
  
  y_double = intt_p19(xx_double);

  % Select only first values
  y = y_double(1:1:N,:);

  % Multiply by 2 to correct for the doubling

  y = mod(y * 2, p);

  
  %% Second method - twist with root of order 2*N to get x^N-1 complex
  
  roots_inv_N2(1,1:2) = [1,0];


  w_inv = [w(1) mod(-w(2),p)];
  
  for i=2:N2
    roots_inv_N2(i,1:2) = power_mod_complex(w_inv, i-1, p);
  end


  y_twisted = intt_p19(xx);

  
  for i=1:N
    y2(i,:) = mult_mod_complex(y_twisted(i,:),roots_inv_N2(i,:),p);
  end

 assert(sum(sum(abs(y-y2)))==0);
  

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
%! y = nega_cyc_ntt_p19(x);
%! z = nega_cyc_intt_p19(y);
%! 
%! assert(sum(sum(abs(z-x)))==0)

