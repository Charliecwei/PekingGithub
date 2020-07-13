function pde=Duichenstokes3data6
%%%3D的stokes方程
%%%u=[10*x*y^4+10*x*z^4-4*x^5,10*y*x^4+10*y*z^4-4*y^5,10*z*x^4+10*z*y^4-4*z^5];

%%% p=-60*x^2*y^2-60*x^2*z^2-60*y^2*z^2+20*(x^4+y^4+z^4);
%%% Dirichlet边界条件


exactDu = struct('u1',@u1,'u2',@u2,'u3',@u3);
pde = struct('f', @f, 'exactp', @exactp, 'exactu',@exactu,'exactDu',exactDu,'g_D',@exactu,'g_N',@g_N);

%%子函数

    function s=f(p)
       x=p(:,1); y=p(:,2);z=p(:,3);
       
       s = [ 160.*x.^3 - 240.*x.*y.^2 - 240.*x.*z.^2, - 240.*x.^2.*y + 160.*y.^3 - 240.*y.*z.^2, - 240.*x.^2.*z - 240.*y.^2.*z + 160.*z.^3];
    end


    function s=exactp(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         
         s=20.*x.^4 - 60.*x.^2.*y.^2 - 60.*x.^2.*z.^2 + 20.*y.^4 - 60.*y.^2.*z.^2 + 20.*z.^4 + 8;
    end


    function s=exactu(p)
         x=p(:,1); y=p(:,2);z=p(:,3);   
         s = [ - 4.*x.^5 + 10.*x.*y.^4 + 10.*x.*z.^4, 10.*x.^4.*y - 4.*y.^5 + 10.*y.*z.^4, 10.*x.^4.*z + 10.*y.^4.*z - 4.*z.^5];
    end


    function s=u1(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         
        s = [ - 20.*x.^4 + 10.*y.^4 + 10.*z.^4,                 40.*x.*y.^3,                 40.*x.*z.^3];
    end



    function s=u2(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         
         s = [                   40.*x.^3.*y, 10.*x.^4 - 20.*y.^4 + 10.*z.^4,                 40.*y.*z.^3];
    end

    function s=u3(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         
         s = [                   40.*x.^3.*z,                 40.*y.^3.*z, 10.*x.^4 + 10.*y.^4 - 20.*z.^4];
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
        
        gg = (abs(-1-x)<1e-14);
        s(gg,1) =-uD1(gg,1)+up(gg);
        s(gg,2) =-uD2(gg,1);
        s(gg,3) =-uD3(gg,1);
        
        gg = (abs(1-x)<1e-14);
        s(gg,1) = uD1(gg,1)-up(gg);
        s(gg,2) = uD2(gg,1);
        s(gg,3) = uD3(gg,1);
        
        
        
        gg = (abs(-1-y)<1e-14);
        s(gg,1) = -uD1(gg,2);
        s(gg,2) = -uD2(gg,2)+up(gg);
        s(gg,3) = -uD3(gg,2);
        
        gg = (abs(1-y)<1e-14);
        s(gg,1) = uD1(gg,2);
        s(gg,2) = uD2(gg,2)-up(gg);
        s(gg,3) = uD3(gg,2);
        
        
        
        
        gg = (abs(-1-z)<1e-14);
        s(gg,1) = -uD1(gg,3);
        s(gg,2) = -uD2(gg,3);
        s(gg,3) = -uD3(gg,3)+up(gg);
        
        gg = (abs(1-z)<1e-14);
        s(gg,1) = uD1(gg,3);
        s(gg,2) = uD2(gg,3);
        s(gg,3) = uD3(gg,3)-up(gg);
        
    end
             


end
         