clear all

po_13 = 13;

r_13 = int64([638, 3892]);

p_13 = int64(2)^po_13-1;

power_mod_complex(r_13,  2^(po_13+1), p_13)

N_pow = 9;

N = 2^N_pow;

r_N_13 = power_mod_complex(r_13, 2^(po_13+1-N_pow), p_13)


  
power_mod_complex(r_N_13,  2^N_pow, p_13)

power_mod_complex(r_N_13,  2^(N_pow-1), p_13)

power_mod_complex(r_N_13,  2^(N_pow-2), p_13)

power_mod_complex(r_N_13,  2^(N_pow-3), p_13)  


ind =  [0   4   2   6   7   3   1   5];

ind_d = 4*ind + 1;

M_13 = int64(zeros(16,16));


om = power_mod_complex(r_N_13,ind_d(1),p_13)

for ll = 1:8
  om = power_mod_complex(r_N_13,ind_d(ll),p_13)
  for k= 1:16
    tmp = power_mod_complex(om, k-1, p_13);
    M_13_real(ll,k) = tmp(1);
    M_13_imag(ll,k) = tmp(2);

    M_13(2*(ll-1)+1,k) = tmp(1);
    M_13(2*(ll-1)+2,k) = tmp(2);

    
    M_inv_13(2*(ll-1)+1,k) = tmp(1);
    M_inv_13(2*(ll-1)+2,k) = p_13-tmp(2);
  end
end




po_19 = 19;

r_19 = int64([502746   91912]);

p_19 = int64(2)^po_19-1;

power_mod_complex(r_19,  2^(po_19+1), p_19)

N_pow = 9;

N = 2^N_pow;

r_N_19 = power_mod_complex(r_19, 2^(po_19+1-N_pow), p_19)


  
power_mod_complex(r_N_19,  2^N_pow, p_19)

power_mod_complex(r_N_19,  2^(N_pow-1), p_19)

power_mod_complex(r_N_19,  2^(N_pow-2), p_19)

power_mod_complex(r_N_19,  2^(N_pow-3), p_19)  


ind =  [0   4   2   6   7   3   1   5];

ind_d = 4*ind + 1;

M_19 = int64(zeros(16,16));

om = power_mod_complex(r_N_19,ind_d(1),p_19)

for ll = 1:8
  om = power_mod_complex(r_N_19,ind_d(ll),p_19)
  for k= 1:16
    tmp = power_mod_complex(om, k-1, p_19);
    M_19_real(ll,k) = tmp(1);
    M_19_imag(ll,k) = tmp(2);

    M_19(2*(ll-1)+1,k) = tmp(1);
    M_19(2*(ll-1)+2,k) = tmp(2);

    M_inv_19(2*(ll-1)+1,k) = tmp(1);
    M_inv_19(2*(ll-1)+2,k) = p_19-tmp(2);
   
  end
end


M = M_13 + 2^13*M_19


M_inv = M_inv_13 + 2^13*M_inv_19


%%csvwrite("M_roots_512.csv",M)
%%csvwrite("M_inv_roots_512.csv",M_inv)


csvwrite("m.csv",reshape(transpose(M),1,256));
%%csvwrite("M_inv_roots_512.csv",reshape(transpose(M_inv),1,256));
