function  [S,v_star_idx] = MAWS_h_theta(u,h,G)
%% the . MA^{WS}_{h,theta} on points (ih,jh)
%% v_star_idx store the v_star index in G

L = size(G,1); 

if nargout == 1 %%Euler method
    for l = 1:L
            e = G(l,:);
            Delta_eu = Delta_e(u,e,h);

            e_vertical = [-e(2) e(1)];    
            Delta_eu_vertical = Delta_e(u,e_vertical,h);
        if l==1
            S = (U_plus(Delta_eu).*U_plus( Delta_eu_vertical));
        else
            S = min([S;(U_plus(Delta_eu).*U_plus( Delta_eu_vertical))]);
        end

    end
elseif nargout ==2 %%Netwon method and find v^*
    S = zeros(L,size(u,1));
    for l = 1:L
        e = G(l,:);
        Delta_eu = Delta_e(u,e,h);

        e_vertical = [-e(2) e(1)];    
        Delta_eu_vertical = Delta_e(u,e_vertical,h);
        
        S(l,:) = (U_plus(Delta_eu).*U_plus( Delta_eu_vertical));
    end
    [S,v_star_idx] = min(S);
    
    v_star_idx = v_star_idx';
end

S = S';
