function pde=Duichenstokes3data1
%%%3D的stokes方程 u=[-y^2(1-y^2),-z^2(1-z^2),-x^2(1-x^2)];
%%% p=0;
%%% Dirichlet边界条件


exactDu = struct('u1',@u1,'u2',@u2,'u3',@u3);
pde = struct('f', @f, 'exactp', @exactp, 'exactu',@exactu,'exactDu',exactDu,'g_D',@exactu,'g_N',@g_N);

%%子函数

    function s=f(p)
       x=p(:,1); y=p(:,2);z=p(:,3);
       
       s(:,1)=2 - 12.*y.^2;
       s(:,2)=2 - 12.*z.^2;
       s(:,3)=2 - 12.*x.^2;
       
    end


    function s=exactp(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         
         s=0.*(x+y+z);
    end


    function s=exactu(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         

        
        s(:,1) = y.^4-y.^2;
        s(:,2) = z.^4-z.^2;
        s(:,3) = x.^4-x.^2;
    end


    function s=u1(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         
         s(:,1)=0*x;
         s(:,2)=(4*y.^3-2*y);
         s(:,3)=0*z;
    end



    function s=u2(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         
         s(:,1)=0*x;
         s(:,2)=0*y;
         s(:,3)=(4*z.^3-2*z);
    end

    function s=u3(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         
         s(:,1)=(4*x.^3-2*x);
         s(:,2)=0*y;
         s(:,3)=0*z;
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
         