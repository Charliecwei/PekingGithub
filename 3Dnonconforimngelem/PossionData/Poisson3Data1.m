function pde=Poisson3Data1
%%%3D的stokes方程 u= z^5*exp(3*x)*cos(y);
%%% Dirichlet边界条件



pde = struct('f', @f,  'exactu',@exactu,'exactDu',@exactDu,'g_D',@exactu,'g_N',@g_N);

%%子函数

    function s=f(p)
       x=p(:,1); y=p(:,2);z=p(:,3);
       s = - 8.*exp(3*x).*cos(y).*z.^5 - 20.*exp(3.*x).*cos(y).*z.^3;     
    end





    function s=exactu(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         s = exp(3*x).*cos(y).*z.^5;
    end




    function s = exactDu(p)
        x = p(:,1); y = p(:,2); z = p(:,3);
        
        s(:,1) = 3*z.^5.*exp(3*x).*cos(y);
        s(:,2) = -z.^5.*exp(3*x).*sin(y);
        s(:,3) = 5*z.^4.*exp(3*x).*cos(y);
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
         