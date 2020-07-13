function phi = getphi2DP2N(lambda)
%%P2非协调元基函数
    phi(:,6) = 4*lambda(:,1).*lambda(:,2);
    phi(:,1) = lambda(:,1).*(2*lambda(:,1)-1);
    phi(:,2) = lambda(:,2).*(2*lambda(:,2)-1);
    phi(:,3) = lambda(:,3).*(2*lambda(:,3)-1);
    phi(:,4) = 4*lambda(:,2).*lambda(:,3);
    phi(:,5) = 4*lambda(:,3).*lambda(:,1);
    phi(:,7) = 2-3*(lambda(:,1).^2+lambda(:,2).^2+lambda(:,3).^2);
end