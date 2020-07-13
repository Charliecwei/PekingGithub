function u = Eulerexplicit(u,h,G)
global f Ftype;
N = 1/h;
L = size(G,1);
for l = 1:L
       e = G(l,:);
        Delta_eu = Delta_e(u,e,h);

        e_vertical = [-e(2) e(1)];    
        Delta_eu_vertical = Delta_e(u,e_vertical,h);
        
        
    if l==1
      
        S = max(U_plus(Delta_eu)+U_plus(Delta_eu_vertical));
        
    else
        
        S = max([S,max(U_plus(Delta_eu)+U_plus(Delta_eu_vertical))]);
                
        
    end
end
S = 1/S;
dt = 0.5*h^2*S;
[i,j] = GlobaltoLocalidx((1:length(u))',N);

switch Ftype
    case 'MAWS'
            Fu = MAWS_h_theta(u,h,G);
    case 'MAWS_delta'
            Fu = MAWS_h_theta_delta(u,h,G);
end

u = u+dt*(Fu-f(i*h,j*h));


end