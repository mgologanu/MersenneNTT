function z = modi(v,p)
  a = v(1);
  b = v(2);

  z = [mod(a,p), mod(b,p)];
  
  %%  z=mod(real(v),p))+mod(imag(v)),p)*i;
end
