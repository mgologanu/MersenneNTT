

input: int32_t a[256]
       uint32_t m[256]
       uint32_t ntt_b[256]

output: int32_t conv_a_b[256]



  load a in A[16,16] line by line

  load m in M[16,16] line by line
  
  load ntt_b in NTT_B[16,16] line by line 
  
  A1 = real_ntt for each column of A
  
  A2 = multiply_complex A1 .* M
  
  NTT_A = complex_ntt for each pair of lines of A2
  
  AB = multiply_complex NTT_A .* NTT_B

  A3 = complex_intt for each pair of lines of AB

  A4 = multiply_complex A3 .* conjugate(M)

  A5 = real_intt for each column of A4

  save A5[16,16] line by line in conv_a_b[256]
       




