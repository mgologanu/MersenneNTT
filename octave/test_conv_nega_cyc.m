clear all;

v = uint32((1:625).^2);

rand ("state", v);


N = 256;

p = 3329;

a = round((2*rand(1,N)-1)*(p-1)/2);

b = round((2*rand(1,N)-1)*(p-1)/2);

tt = my_conv(a,b);
abz(1:N) = tt(1:N);
abz(1:N-1) = abz(1:N-1)-tt(N+1:2*N-1);


ab1 = conv_mod_nega_cyc(a,b,p);




p13 = int64(2)^13-1;

p19 = int64(2)^19-1;

a13 = mod(a,p13);
b13 = mod(b,p13);


A13 = nega_cyc_ntt_p13(a13);
B13 = nega_cyc_ntt_p13(b13);


AB13 = mult_mod_complex(A13,B13,int64(2)^13-1);


ab13 = nega_cyc_intt_p13(AB13);

ab13 = ab13(:,1);

ab13'


p19 = int64(2)^19-1;

a19 = mod(a,p19);
b19 = mod(b,p19);


A19 = nega_cyc_ntt_p19(a19);
B19 = nega_cyc_ntt_p19(b19);


AB19 = mult_mod_complex(A19,B19,int64(2)^19-1);


ab19 = nega_cyc_intt_p19(AB19);

ab19 = ab19(:,1);

ab19'




zz = int64(ab13') + mod((int64(ab19')-int64(ab13'))*int64(8321),p19)*p13;



ind = zz > (p13*p19-1)/2;

zz(ind) = zz(ind) - p13*p19;

assert(sum(abs(abz - zz))==0)


ab2 = mod(zz,p);

assert(sum(abs(ab1-ab2))==0)

ind1 =  [0   4   2   6   1   5   3   7]

ind2 =  [0   4   2   6   7   3   1   5]

ind_complex = [2*ind1 2*ind2+1];

ind_real =  [0   4   2   6   7   3   1   5];


tmp=transpose(reshape((0:127)',16,8));
tmp1=tmp(ind_real+1,:);
tmp2=tmp1(:,ind_complex+1);

%%ind_B=reshape(transpose(tmp2),1,128)



B13_re = B13(1:128,1);
B13_im = B13(1:128,2);

B13_reordered_re = B13_re(tmp2+1);
B13_reordered_im = B13_im(tmp2+1);

for k = 1:8
  B13_reordered(2*(k-1)+1,:) = B13_reordered_re(k,:);
  B13_reordered(2*(k-1)+2,:) = B13_reordered_im(k,:);
end



B19_re = B19(1:128,1);
B19_im = B19(1:128,2);

B19_reordered_re = B19_re(tmp2+1);
B19_reordered_im = B19_im(tmp2+1);

for k = 1:8
  B19_reordered(2*(k-1)+1,:) = B19_reordered_re(k,:);
  B19_reordered(2*(k-1)+2,:) = B19_reordered_im(k,:);
end


B_reordered = B13_reordered + 2^13 * B19_reordered;



csvwrite("ntt_b.csv",reshape(transpose(B_reordered),1,256));

csvwrite("a.csv",a);


ind_ab = ab2 > (p-1)/2;

ab_signed = ab2;
ab_signed(ind_ab) = ab_signed(ind_ab) - p;


csvwrite("conv_a_b.csv", ab_signed);
