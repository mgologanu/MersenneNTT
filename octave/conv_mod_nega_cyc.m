function  z = conv_mod_nega_cyc(x,y,p);


  N = length(x);
  
  assert(length(y)==N);

  

  
  zz = conv_mod(x,y,p);

  z(1:N) = zz(1:N);

  z(1:N-1) = mod(z(1:N-1) - zz(N+1:2*N-1),p);
  
  
end
