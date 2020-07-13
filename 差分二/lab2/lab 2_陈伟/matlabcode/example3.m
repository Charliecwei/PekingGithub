function pde = example3

pde  = struct('u',@u,'f',@f,'g_D',@u);

function s = f(X)
    x = X(:,1);
    y = X(:,2);
    s = 0*x;
    idx = abs(y)<=abs(x).^3;
    s(idx) = 36-9*y(idx).^2./x(idx).^6;
    s(~idx) = 8/9-5/9*x(~idx).^2./power(y(~idx).^2,1/3);
end

function  s = u(X)
   x = X(:,1);
   y = X(:,2);
   s = 0*x;
   idx = abs(y)<abs(x).^3; 
   s(idx) = x(idx).^4+1.5*y(idx).^2./x(idx).^2;
   s(~idx) = 0.5*x(~idx).^2.*power(y(~idx).^2,1/3)+2*power(y(~idx).^4,1/3);
end

end