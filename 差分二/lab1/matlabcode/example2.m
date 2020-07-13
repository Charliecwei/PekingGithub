function pde = example2

pde  = struct('u',@u,'f',@f,'g_D',@u);

function s = f(x,y)
   s = 1-0.2./sqrt((x-0.5).^2+(y-0.5).^2);
   idx = s<0;
   s(idx) = 0;
end

function  s = u(x,y)
    s = sqrt((x-0.5).^2+(y-0.5).^2)-0.2;
    idx = s<0;
    s(idx) = 0;
    s = 0.5*s.^2;
end

end