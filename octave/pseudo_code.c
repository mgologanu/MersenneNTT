

input: a[256], m[256], ntt_b[256]

output: conv_a_b[256]



  load A[16,16] line by line

  load M[16,16] line by line
  
  load NTT_B[16,16] line by line 
  
  A1 = real_ntt for each column of A
  
  A2 = multiply_complex A1 .* M
  
  NTT_A = complex_ntt for each pair of lines of A2
  
  AB = multiply_complex NTT_A .* NTT_B

  A3 = complex_intt for each pair of lines of AB

  A4 = multiply_complex A3 .* conjugate(M)

  A5 = real_intt for each column of A4

  save A5[16,16] line by line in conv_a_b[256]
       




