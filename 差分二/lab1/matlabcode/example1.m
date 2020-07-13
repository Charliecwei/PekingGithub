function pde = example1

pde  = struct('u',@u,'f',@f,'g_D',@u);

function s = f(x,y)
   s = (1+x.^2+y.^2).*exp(x.^2+y.^2);
end

function  s = u(x,y)
    s = exp((x.^2+y.^2)/2);
end

end