function res = mat_vect_mod_complex(M,x,p)

  %% mat_vect_mod_complex(M,x,p) - returns M * x mod p
  %%
  %% where M is a matrix of complex remainders modulo p and
  %% x is a vector of complex remainders modulo p
  %%
  %% M(:,:,1) is the real part, M(:,:,2) is the imaginary part
  %%
  %% x(:,1) is the real part, x(:,1) is the imaginary part
  
  [rows, cols, plane_M] = size(M);

  [rows_x, plane_x] = size(x);

  if (plane_x ~= plane_M || rows_x ~= cols )
    disp("Wrong input: [rows, cols, plane_M, rows_x, plane_x] ");
    disp([rows, cols, plane_M, rows_x, plane_x]);
    res =  NaN;
    return;
  endif
	   
  for i = 1:rows
    tmp_re = int64(0);
    tmp_im = int64(0);
    for j = 1:cols
      tmp_re = mod(tmp_re + mod(mod(M(i,j,1)*x(j,1),p) - mod(M(i,j,2)*x(j,2),p),p),p) ;
      tmp_im = mod(tmp_im + mod(mod(M(i,j,1)*x(j,2),p) + mod(M(i,j,2)*x(j,1),p),p),p) ;
    end
    res(i,1) = tmp_re;
    res(i,2) = tmp_im;
  end
  

end


%!test
%!
%!
%!
%! p = int64(2)^13-1;
%!
%! re = int64([2600, 2980, 1480, 3429, 6608, 3328, 5255, 5870, 2168, 454, 2764, 1978, 826, 551, 3783, 7693]);
%!
%! im = int64([7661, 3262, 3695, 1096, 8104, 311, 1626, 4582, 2052, 1645, 6541, 5549, 7340, 4572, 7779, 1445]);
%!
%!
%! x(:,1) = re';
%! x(:,2) = im';
%!
%! N = 16;
%!
%! w_16 = int64([6456 812]);
%!
%! roots_16(1,1:2) = [1,0];
%!
%! for i=2:16
%!   roots_16(i,1:2) = power_mod_complex(w_16, i-1, p);
%! end
%! 
%! M=int64(zeros(N,N,2));
%! 
%! kk = int64([1:N]-1)';
%! 
%! for j=1:N
%!   po = mod((j-1)*kk,N) + 1;
%!   M(:,j,1) =roots_16(po,1);
%!   M(:,j,2) =roots_16(po,2);
%! end
%! 
%! MM=int64(zeros(N,N,2));
%! 
%! kk = int64([1:N]-1);
%! 
%! for j=1:N
%!   po = mod((j-1)*kk,N) + 1;
%!   MM(j,:,1) =roots_16(po,1)';
%!   MM(j,:,2) =roots_16(po,2)';
%! end

%! assert(sum(sum(sum(abs(M(:,:,:)-MM(:,:,:)))))==0)

%! y = mat_vect_mod_complex(M,x,p);

%! re_expected = [2621   266  6605  3181  1802   456  4179  7252  7392  3647  2783  4898  4229  2494   342  5835];

%! im_expected = [1732  6168  7129  1754  2050  3203  2648  4458  5954  5278  4585    27   791  6342  3669  1260];

%! y_expected(:,1) = re_expected';
%! y_expected(:,2) = im_expected';

%! assert(sum(sum(abs(y-y_expected)))==0)
