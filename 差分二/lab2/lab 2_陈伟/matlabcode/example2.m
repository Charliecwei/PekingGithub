function pde = example2

pde  = struct('u',@u,'f',@f,'g_D',@u);

function s = f(X)
    x = X(:,1);
    y = X(:,2);
    R = sqrt(x.^2+y.^2);
    idx = (R>=1/2);
    s = 16 + 0*R;
    s(idx) = 64 - 16./R(idx);
    
end

function  s = u(X)
    x = X(:,1);
    y = X(:,2);
    R2 = x.^2+y.^2;
    s = 2*R2;
    idx =( R2>=1/4);
    s(idx) = s(idx)+ 2*power(sqrt(R2(idx))-1/2,2);
    
end

end