function  z = my_conv(x,y);

  ll = length(x) + length(y) - 1;

  for i = 1:ll
    z(i) = 0;
    
    for j = 1:length(x)
      if ((i - j >= 0) && (i - j < length(y))) 
        z(i) = z(i) +  x(j) * y(i - j + 1);
      endif
    endfor
  endfor
  
end
