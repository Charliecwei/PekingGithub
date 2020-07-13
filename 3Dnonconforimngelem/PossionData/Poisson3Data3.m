function pde=Poisson3Data3
%%%3D的stokes方程 u=  (x-1/2).^5.*y^5.*(z-1/2);
%%% Dirichlet边界条件



pde = struct('f', @f,  'exactu',@exactu,'exactDu',@exactDu,'g_D',@exactu,'g_N',@g_N);

%%子函数

    function s=f(p)
       x=p(:,1); y=p(:,2);z=p(:,3);
       s = - 20*y.^3.*(x - 1/2).^5.*(z - 1/2) - 20*y.^5.*(x - 1/2).^3.*(z - 1/2);
    end





    function s=exactu(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         s =  (x-1/2).^5.*y.^5.*(z-1/2);
    end




    function s = exactDu(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        
        s(:,1) = 5*y.^5.*(x - 1/2).^4.*(z - 1/2);
        s(:,2) = 5*y.^4.*(x - 1/2).^5.*(z - 1/2);
        s(:,3) = y.^5.*(x - 1/2).^5;
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
         