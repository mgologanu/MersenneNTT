function res = power_modi(t, n, p)

  %% Calculates t^n mod p

  max_power = floor(log2(n));

  pow2 = int64(zeros(max_power+1,2));

  pow2(1,:) = t;
  
  for k = 1:max_power
    pow2(k+1,:) = modi( mult_c(pow2(k,:),pow2(k,:)), p);
  end

  %%disp(pow2)
  ss=(dec2bin(n));


  res = [double(1), 0];
  
  for k = 1:length(ss)
    if (ss(k) == '1')
     %% disp(k)
     %% disp(pow2(max_power+2-k))
      res = modi(mult_c(res, pow2(max_power+2-k,:)),p);
    end
  end
    
  
end
