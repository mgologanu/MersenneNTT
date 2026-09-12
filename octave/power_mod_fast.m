function res = power_mod_fast(t, n, p)

  %% Calculates t^n mod p

  max_power = floor(log2(n));

  pow2 = int64(zeros(1,max_power+1));

  pow2(1) = int64(t);
  
  for k = 1:max_power
    pow2(k+1) = mod( pow2(k)*pow2(k), p);
  end

  %%disp(pow2);
  
  ss=(dec2bin(n));

  

  res = int64(1);
  
  for k = 1:length(ss)
    if (ss(k) == '1')
     %% disp(k)
     %% disp(pow2(max_power+2-k))
      res = mod(res * pow2(max_power+2-k),p);
    end
  end
    
  
end
