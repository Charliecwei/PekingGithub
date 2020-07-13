function mu = DYNVIS(T)
    %Dynamic viscosity
    global mu_0 T_0;
    mu = mu_0.*power(T./T_0,1.5).*(T_0+110)./(T+110);
end
