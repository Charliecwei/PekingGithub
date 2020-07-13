function pde=Duichenstokes3data3
%%%3D的stokes方程 u=[-y^2(1-y^2),-z^2(1-z^2),-x^2(1-x^2)];
%%% p=0;
%%% Dirichlet边界条件


exactDu = struct('u1',@u1,'u2',@u2,'u3',@u3);
pde = struct('f', @f, 'exactp', @exactp, 'exactu',@exactu,'exactDu',exactDu,'g_D',@exactu,'g_N',@g_N);

%%子函数

    function s=f(p)
       x=p(:,1); y=p(:,2);z=p(:,3);
       
       s(:,1)=6.*y.*z.*exp(x) - (y - 1/2).*(z - 1) - 6.*cos(pi.*z) + 3.*x.^2.*pi^2.*cos(pi.*z) + 6;
       s(:,2)=x.^3.*y.*pi^4.*cos(pi.*z) - 6.*z.*exp(x) - 3.*y.^2.*z.*exp(x) - 12.*x.*y.*pi^2.*cos(pi.*z) - (x - 1/2).*(z - 1);
       s(:,3)=6.*x.*pi.*sin(pi.*z) - y.^3.*exp(x) - 6.*y.*exp(x) - x.^3.*pi^3.*sin(pi.*z) - (x - 1/2).*(y - 1/2);
 
 
       
    end


    function s=exactp(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         
         s=-(x - 1/2).*(y - 1/2).*(z - 1);
    end


    function s=exactu(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         

        
        s(:,1) = 3.*x.^2.*cos(pi.*z) - 3.*(x - 1).^2 - 6.*y.*z.*exp(x);
        s(:,2) = 3.*y.^2.*z.*exp(x) - 6.*x.*y.*cos(pi.*z) + x.^3.*y.*pi.^2.*cos(pi.*z);
        s(:,3) = 3.*z.*(2.*x - 2) + y.^3.*exp(x) - x.^3.*pi.*sin(pi.*z);
    end


    function s=u1(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         
         s = [  6.*x.*cos(pi.*z) - 6.*x - 6.*y.*z.*exp(x) + 6,   -6.*z.*exp(x), - 6.*y.*exp(x) - 3.*x.^2.*pi.*sin(pi.*z)];
    end



    function s=u2(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         
        s = [ 3.*y.^2.*z.*exp(x) - 6.*y.*cos(pi.*z) + 3.*x.^2.*y.*pi.^2.*cos(pi.*z), 6.*y.*z.*exp(x) - 6.*x.*cos(pi.*z) + x.^3.*pi.^2.*cos(pi.*z), 3.*y.^2.*exp(x) + 6.*x.*y.*pi.*sin(pi.*z) - x.^3.*y.*pi.^3.*sin(pi.*z)];

 
    end

    function s=u3(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         
         s = [    6.*z + y.^3.*exp(x) - 3.*x.^2.*pi.*sin(pi.*z),  3.*y.^2.*exp(x), - pi.^2.*cos(z.*pi).*x.^3 + 6.*x - 6];
    end



  function s = g_N(p)
        x=p(:,1); y=p(:,2);z=p(:,3);
        
       
        uD(:,:,1) = u1(p);
        uD(:,:,2) = u2(p);
        uD(:,:,3) = u3(p);
        
        uD1 = zeros(length(x),3);
        uD2 = zeros(length(x),3);
        uD3 = zeros(length(x),3);
        
        for i = 1:3
            uD1(:,i) = uD(:,i,1)+uD(:,1,i);
            uD2(:,i) = uD(:,i,2)+uD(:,2,i);
            uD3(:,i) = uD(:,i,3)+uD(:,3,i);
        end
        
        
        up = exactp(p);
        
        s = zeros( size(p,1),3);
        
        gg = (abs(x)<=1e-15);
        s(gg,1) =-uD1(gg,1)+up(gg);
        s(gg,2) =-uD2(gg,1);
        s(gg,3) =-uD3(gg,1);
        
        gg = (abs(1-x)<=1e-15);
        s(gg,1) = uD1(gg,1)-up(gg);
        s(gg,2) = uD2(gg,1);
        s(gg,3) = uD3(gg,1);
        
        
        
        gg = (abs(y)<=1e-15);
        s(gg,1) = -uD1(gg,2);
        s(gg,2) = -uD2(gg,2)+up(gg);
        s(gg,3) = -uD3(gg,2);
        
        gg = (abs(1-y)<=1e-15);
        s(gg,1) = uD1(gg,2);
        s(gg,2) = uD2(gg,2)-up(gg);
        s(gg,3) = uD3(gg,2);
        
        
        
        
        gg = (abs(z)<=1e-15);
        s(gg,1) = -uD1(gg,3);
        s(gg,2) = -uD2(gg,3);
        s(gg,3) = -uD3(gg,3)+up(gg);
        
        gg = (abs(1-z)<=1e-15);
        s(gg,1) = uD1(gg,3);
        s(gg,2) = uD2(gg,3);
        s(gg,3) = uD3(gg,3)-up(gg);
        
    end
        


end
         