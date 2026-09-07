function y = nega_cyc_ntt_p13(x);

  
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

  
  if (2^N_pow ~= N || N_pow > 13)
    disp("Length of array x should be a power of 2 no larger than 2^13");
    y = NaN;
    return;
  endif


  N2 = 2*N;

  N2_pow = N_pow + 1;
  
  p_pow = 13;
  
  p = int64(2)^p_pow-1;

  w0 = int64([638, 3892]);

  w = power_mod_complex(w0, 2^(p_pow+1 - N2_pow), p);
 
  sq2 = 2^((p_pow-1)/2);
  
  assert(sum(abs([  1   0] - power_mod_complex(w,2^(N2_pow),  p)))==0)
  assert(sum(abs([p-1   0] - power_mod_complex(w,2^(N2_pow-1),p)))==0)
  assert(sum(abs([  0   1] - power_mod_complex(w,2^(N2_pow-2),p)))==0)
  assert(sum(abs([sq2 sq2] - power_mod_complex(w,2^(N2_pow-3),p)))==0)
  
  

  %% First method - double the original array [x -x]
  
  xx_double = xx;
  
  xx_double(N+1:2*N, :) = mod(-xx_double(1:N, :),p);


  y_double = ntt_p13(xx_double);

  % Select only non-zero values
  y = y_double(2:2:N2,:);

  % Divide by 2 to correct for the doubling
  two_inv = int64(2)^(p_pow-1);
    
  y = mod(y * two_inv, p);

  
  %% Second method - twist coefficients with root of order 2*N to get
  %% x^N-1 complex
  
  roots_N2(1,1:2) = [1,0];

  for i=2:N2
    roots_N2(i,1:2) = power_mod_complex(w, i-1, p);
  end
  
  for i=1:N
    xx_twisted(i,:) = mult_mod_complex(xx(i,:),roots_N2(i,:),p);
  end

  y_twisted = ntt_p13(xx_twisted);


  assert(sum(sum(abs(y-y_twisted)))==0);
  

end

  
