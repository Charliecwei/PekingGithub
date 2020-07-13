function Fu = MAWS_h_theta_delta(u,h,G)
L = size(G,1); 

for l = 1:L
        e = G(l,:);
        Delta_eu = Delta_e(u,e,h);

        e_vertical = [-e(2) e(1)];    
        Delta_eu_vertical = Delta_e(u,e_vertical,h);
    if l==1
        Fu = ( U_plus_delta(Delta_eu).*U_plus_delta( Delta_eu_vertical));
    else
        Fu = min_delta(Fu,(U_plus_delta(Delta_eu).*U_plus_delta( Delta_eu_vertical)));
    end

end