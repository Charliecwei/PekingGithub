function pde = example1

pde  = struct('u',@u,'f',@f,'g_D',@u);

function s = f(X)
    x = X(:,1);
    y = X(:,2);
   s = (1+x.^2+y.^2).*exp(x.^2+y.^2);
end

function  s = u(X)
    x = X(:,1);
    y = X(:,2);
    s = exp((x.^2+y.^2)/2);
end

end