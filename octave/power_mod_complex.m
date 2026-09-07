function res = power_mod_complex(t, n, p)

  %% power_mod_complex(t, n, p) - returns t^n mod p
  %%
  %% where t = [real_part imag_part] is a complex remainder modulo p
  

  max_power = floor(log2(n));

  pow2 = int64(zeros(max_power+1,2));

  pow2(1,:) = t;
  
  for k = 1:max_power
    pow2(k+1,:) = mult_mod_complex(pow2(k,:),pow2(k,:), p);
  end

  ss=(dec2bin(n));


  res = int64([1, 0]);
  
  for k = 1:length(ss)
    if (ss(k) == '1')
      res = mult_mod_complex(res, pow2(max_power+2-k,:),p);
    end
  end
    
  
end

%!test
%!
%! p = int64(2)^13-1;
%! r = int64([638, 3892]);
%! assert (norm (double([1 0] - power_mod_complex(r,2^14,p)))==0)
%! assert (norm (double([p-1 0] - power_mod_complex(r,2^13,p)))==0)
%! assert (norm (double([0 1] - power_mod_complex(r,2^12,p)))==0)
%! assert (norm (double([64 64] - power_mod_complex(r,2^11,p)))==0)



%!test
%!
%! p = int64(2)^19-1;
%! r = int64([502746   91912]);
%! assert (norm (double([1 0] - power_mod_complex(r,2^20,p)))==0)
%! assert (norm (double([p-1 0] - power_mod_complex(r,2^19,p)))==0)
%! assert (norm (double([0 1] - power_mod_complex(r,2^18,p)))==0)
%! assert (norm (double([512 512] - power_mod_complex(r,2^17,p)))==0)
