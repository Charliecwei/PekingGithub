function [Fu,Grad_Fu] = Grad_MAWS_h_theta(u,h,G)
%%To calculate the Grad_u of MAWS_h_theta and Fu
Nu = length(u);
%N = 1/h;

[Fu,v_star_idx]=MAWS_h_theta(u,h,G);

 v_star = G(v_star_idx,:);
 
 v_star_vertical = [-v_star(:,2),v_star(:,1)];
 


 ii = (1:Nu)';
 jj = ii;
 
 Delta_eu = Delta_e(u,v_star,h);
 
 Delta_eu_v = Delta_e(u,v_star_vertical,h);
 
 GraduDeltave = Grad_uOfDelta_eu(u,v_star,h);
 
 GraduDeltave_v = Grad_uOfDelta_eu(u,v_star_vertical,h);
 
 Grad_Fu =  sparse(ii,jj,Delta_eu,Nu,Nu)*GraduDeltave_v...
                    + sparse(ii,jj,Delta_eu_v,Nu,Nu)*GraduDeltave;

               
 

    
end
    
