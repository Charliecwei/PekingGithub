function pde=Poisson3Data5
%%%3D的Poisson方程 u=  sin(2*pi*x)*sin(2*pi*y)*sin(2*pi*z);
%%% Dirichlet边界条件



pde = struct('f', @f,  'exactu',@exactu,'exactDu',@exactDu,'g_D',@exactu,'g_N',@g_N);

%%子函数

    function s=f(p)
       x=p(:,1); y=p(:,2);z=p(:,3);
       s = 12.*pi.^2.*sin(2*pi*x).*sin(2*pi*y).*sin(2*pi*z);
    end





    function s=exactu(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         s =  sin(2*pi*x).*sin(2*pi*y).*sin(2*pi*z);
    end




    function s = exactDu(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        
        s(:,1) = 2*pi*cos(2*pi*x).*sin(2*pi*y).*sin(2*pi*z);
        s(:,2) = 2*pi*cos(2*pi*y).*sin(2*pi*x).*sin(2*pi*z);
        s(:,3) = 2*pi*cos(2*pi*z).*sin(2*pi*x).*sin(2*pi*y);
    end



  function s = g_N(p)
        x=p(:,1); y=p(:,2);z=p(:,3);
        
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
       
        
        
        
        
        gg = (abs(z)<=1e-15);
        s(gg) = -uD(gg,3);
   
        
        gg = (abs(1-z)<=1e-15);
        s(gg) = uD(gg,3);
        
    end
        


end
         