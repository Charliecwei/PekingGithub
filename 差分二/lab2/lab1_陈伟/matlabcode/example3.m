function pde = example3

pde  = struct('u',@u,'f',@f,'g_D',@u);

function s = f(x,y)
   s = 2./((2-(x.^2+y.^2)).^2);
end

function  s = u(x,y)
    s = -sqrt(2-(x.^2+y.^2));
end

end