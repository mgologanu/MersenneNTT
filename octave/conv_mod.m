function  z = conv_mod(x,y,p);

  ll = length(x) + length(y) - 1;

  for i = 1:ll
    z(i) = 0;
    
    for j = 1:length(x)
      if ((i - j >= 0) && (i - j < length(y))) 
        z(i) = mod(z(i) +  mod(x(j) * y(i - j + 1),p),p);
      endif
    endfor
  endfor
  
end
