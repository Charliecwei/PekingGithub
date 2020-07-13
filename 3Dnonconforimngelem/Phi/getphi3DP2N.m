function phi = getphi3DP2N(lambda)
%%P2非协调元基函数in 3D

    phi(:,10)= 4*lambda(:,3).*lambda(:,4);
    phi(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
    phi(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
    phi(:,3) = lambda(:,3).*(2*lambda(:,3)-1);
    phi(:,4) = lambda(:,4).*(2*lambda(:,4)-1);
    phi(:,5) = 4*lambda(:,1).*lambda(:,2);
    phi(:,6) = 4*lambda(:,1).*lambda(:,3);
    phi(:,7) = 4*lambda(:,1).*lambda(:,4);
    phi(:,8) = 4*lambda(:,2).*lambda(:,3);
    phi(:,9) = 4*lambda(:,2).*lambda(:,4);
    phi(:,11) = 2-4*(lambda(:,1).^2+lambda(:,2).^2+lambda(:,3).^2+lambda(:,4).^2);
   % phi(:,11) = lambda(:,1).*lambda(:,2).*lambda(:,3).*lambda(:,4);
end