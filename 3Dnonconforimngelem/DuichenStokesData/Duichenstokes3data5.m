function pde=Duichenstokes3data5
%%%3D��stokes���� u=[ pi*cos(pi*y)*sin(pi*x)^2*sin(pi*y)*sin(pi*z),
%%%-pi*cos(pi*x)*sin(pi*x)*sin(pi*y)^2*sin(pi*z), 0];

%%% p=0;
%%% Dirichlet�߽�����


exactDu = struct('u1',@u1,'u2',@u2,'u3',@u3);
pde = struct('f', @f, 'exactp', @exactp, 'exactu',@exactu,'exactDu',exactDu,'g_D',@exactu,'g_N',@g_N);

%%�Ӻ���

    function s=f(p)
       x=p(:,1); y=p(:,2);z=p(:,3);
       
       s = [ 7.*pi.^3.*cos(pi.*y).*sin(pi.*x).^2.*sin(pi.*y).*sin(pi.*z) - 2.*pi.^3.*cos(pi.*x).^2.*cos(pi.*y).*sin(pi.*y).*sin(pi.*z), 2.*pi.^3.*cos(pi.*x).*cos(pi.*y).^2.*sin(pi.*x).*sin(pi.*z) - 7.*pi.^3.*cos(pi.*x).*sin(pi.*x).*sin(pi.*y).^2.*sin(pi.*z), 0*x];

    end


    function s=exactp(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         
         s=0*(x+y+z);
    end


    function s=exactu(p)
         x=p(:,1); y=p(:,2);z=p(:,3);   
         s = [ pi.*cos(pi.*y).*sin(pi.*x).^2.*sin(pi.*y).*sin(pi.*z), -pi.*cos(pi.*x).*sin(pi.*x).*sin(pi.*y).^2.*sin(pi.*z), 0*z];
    end


    function s=u1(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         
        s = [                        2.*pi.^2.*cos(pi.*x).*cos(pi.*y).*sin(pi.*x).*sin(pi.*y).*sin(pi.*z), pi.^2.*cos(pi.*y).^2.*sin(pi.*x).^2.*sin(pi.*z) - pi.^2.*sin(pi.*x).^2.*sin(pi.*y).^2.*sin(pi.*z),  pi.^2.*cos(pi.*y).*cos(pi.*z).*sin(pi.*x).^2.*sin(pi.*y)];
    end



    function s=u2(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         
         s = [ pi.^2.*sin(pi.*x).^2.*sin(pi.*y).^2.*sin(pi.*z) - pi.^2.*cos(pi.*x).^2.*sin(pi.*y).^2.*sin(pi.*z),                       -2.*pi.^2.*cos(pi.*x).*cos(pi.*y).*sin(pi.*x).*sin(pi.*y).*sin(pi.*z), -pi.^2.*cos(pi.*x).*cos(pi.*z).*sin(pi.*x).*sin(pi.*y).^2];
    end

    function s=u3(p)
         x=p(:,1); y=p(:,2);z=p(:,3);
         
         s = [0*x,0*y,0*z];
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
        
        gg = (abs(x)<1e-14);
        s(gg,1) =-uD1(gg,1)+up(gg);
        s(gg,2) =-uD2(gg,1);
        s(gg,3) =-uD3(gg,1);
        
        gg = (abs(1-x)<1e-14);
        s(gg,1) = uD1(gg,1)-up(gg);
        s(gg,2) = uD2(gg,1);
        s(gg,3) = uD3(gg,1);
        
        
        
        gg = (abs(y)<1e-14);
        s(gg,1) = -uD1(gg,2);
        s(gg,2) = -uD2(gg,2)+up(gg);
        s(gg,3) = -uD3(gg,2);
        
        gg = (abs(1-y)<1e-14);
        s(gg,1) = uD1(gg,2);
        s(gg,2) = uD2(gg,2)-up(gg);
        s(gg,3) = uD3(gg,2);
        
        
        
        
        gg = (abs(z)<1e-14);
        s(gg,1) = -uD1(gg,3);
        s(gg,2) = -uD2(gg,3);
        s(gg,3) = -uD3(gg,3)+up(gg);
        
        gg = (abs(1-z)<1e-14);
        s(gg,1) = uD1(gg,3);
        s(gg,2) = uD2(gg,3);
        s(gg,3) = uD3(gg,3)-up(gg);
        
    end
             


end
         