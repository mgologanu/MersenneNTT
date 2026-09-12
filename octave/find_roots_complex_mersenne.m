clear all;

ex = 19 

p = int64(2^ex-1)

q = int64(2^(ex-2))

factors = unique(sort(factor(p-1)))


for tt = 7:100
    ok = true;
    for k = 1:length(factors)
        temp = power_mod_fast(tt,(p-1)/factors(k),p)
        if (temp==1)
            ok = false;
            break;
        end
    end
    if(ok)
        break
    end
end


assert(ok)
primitive = int64(tt)




a = power_mod_fast(int64(2),q,p);
b = power_mod_fast(int64(3),q,p);

x = [a,b]

for k = 1:ex+1
    x = mult_mod_complex(x,x,p);
end


for l =2:1:p-1
    x = int64([1,l]);
    for k = 1:ex
        x = mult_mod_complex(x,x,p);
    end
    if (x(2)==0 && x(1)==1+l*l)

        no = 1+l^2

        primitive_sq = mod(primitive*primitive,p);

        y = primitive;
        ok = false;
        for ll=1:(p-1)/2
            if (y == no)
            	ok = true;
            	break;
            end
            y = mod(y*primitive_sq,p);
        end

        if(ok)
            break
        end
    end
end


po = 2*ll-1

if(true)

    primitive_inv = power_mod_fast(primitive,p-2,p)

    assert(mod(primitive_inv*primitive,p)==1);

    temp = power_mod_fast(primitive_inv,(po-1)/2,p);

    c = modi([temp,l*temp],p)

    po_target = 2^(ex+1)

    x=power_mod_complex(c,(p^2-1)/po_target,p)


end

