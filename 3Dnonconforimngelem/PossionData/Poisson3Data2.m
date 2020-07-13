function pde=Poisson3Data2
%%%3D的stokes方程 u= -y^10*z^9*(x - 1)^5;
%%% Dirichlet边界条件



pde = struct('f', @f,  'exactu',@exactu,'exactDu',@exactDu,'g_D',@exactu,'g_N',@g_N);

%%子函数

    function s=f(p)
       x=p(:,1); y=p(:,2);z=p(:,3);
       s = 90.*y.^8.*z.^9.*(x - 1).^5 + 72.*y.^10.*z.^7.*(x - 1).^5 + 20*y.^10.*z.^9.*(x - 1).^3;
    end





    function s=exactu(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         s = -y.^10.*z.^9.*(x - 1).^5;
    end




    function s = exactDu(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        
        s(:,1) = -5.*y.^10.*z.^9.*(x - 1).^4;
        s(:,2) = -10.*y.^9.*z.^9.*(x - 1).^5;
        s(:,3) = -9.*y.^10.*z.^8.*(x - 1).^5;
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
         