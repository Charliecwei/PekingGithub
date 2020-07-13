function pde=Poisson2DData2
%%%3D的stokes方程 u= exp(5*x+y)
%%% Dirichlet边界条件



pde = struct('f', @f,  'exactu',@exactu,'exactDu',@exactDu,'g_D',@exactu,'g_N',@g_N);

%%子函数

    function s=f(p)
       x=p(:,1); y=p(:,2);
       s = 2*pi^2*sin(pi*x).*sin(pi*y);
    end





    function s=exactu(p)
         x=p(:,1); y=p(:,2);
         s = sin(pi*x).*sin(pi*y);
    end




    function s = exactDu(p)
        x = p(:,1); y = p(:,2);
        
        s(:,1) = pi*cos(pi*x).*sin(pi*y);
        s(:,2) = pi*sin(pi*x).*cos(pi*y);

    end



  function s = g_N(p)
        x=p(:,1); y=p(:,2);
        
        uD = exactDu(p);
        
        s = zeros( size(p,1),1);
        
        gg = (abs(x)<=1e-15);
        s(gg) =-uD(gg,1);
        
        gg = (abs(1-x)<=1e-15);
        s(gg) = uD(gg,1);
        
        
        
        gg = (abs(y)<=1e-15);
        s(gg) = -uD(gg,2);
        
        gg = (abs(1-y)<=1e-15);
        s(gg) = uD(gg,2);
       
        
        
        
    end
        


end
         