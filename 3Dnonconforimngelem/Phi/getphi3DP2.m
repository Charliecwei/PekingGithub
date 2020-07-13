function phi = getphi3DP2(lambda)
%%P2Ôª»ùº¯Êý
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
end