function res = mult_c(x,y)
  a0 = x(1);
  b0 = x(2);

  a1 = y(1);
  b1 = y(2);

  res = [a0*a1 - b0*b1, a0*b1 + a1*b0];


end
