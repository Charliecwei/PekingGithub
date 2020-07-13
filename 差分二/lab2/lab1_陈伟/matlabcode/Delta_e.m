function [Delta_eu, Gradu_Deltaeu] = Delta_e(u,e,h)
%%Find the Delta of u at (ih,jh) in the direction of e, and the grid size is h
global g_D;
N = 1/h;


u_l = zeros(size(u));
u_r = zeros(size(u));

rho_l = ones(size(u));
rho_r = ones(size(u));

[i,j] = GlobaltoLocalidx((1:length(u))',N);

if size(e,1) ==1
    e = repmat(e,size(u));
end





    i_l = i-rho_l.*e(:,1);%the x value of the interpolation point on the left
    j_l = j-rho_l.*e(:,2);%the y value of the interpolation point on the left




    i_r = i+rho_r.*e(:,1);%the x value of the interpolation point on the left
    j_r = j+rho_r.*e(:,2);%the y value of the interpolation point on the left


%% left

    % find the inner point
    idx = and(and(1<=i_l,i_l<=N-1),and(1<=j_l,j_l<=N-1));
    s = LocaltoGlobalidx(i_l(idx),j_l(idx),N);
    u_l(idx) = u(s);  
    
    
    %boundary point
    idx_bd = ~idx;
    e_bd = e(idx_bd,:);
    
    i_bd = i(idx_bd);
    j_bd = j(idx_bd);
    rho_l_bd = rho_l(idx_bd);
   
    
    
    
    
    
    rho_i_r = (i_bd-N)./e_bd(:,1);
    idx = rho_i_r<0;
    rho_i_r(idx) = zeros(size(rho_i_r(idx)))+1;
    
    rho_i_l = (i_bd-0)./e_bd(:,1);
    idx = rho_i_l<0;
    rho_i_l(idx) = zeros(size(rho_i_l(idx)))+1;
    
    
   
    rho_j_r = (j_bd-N)./e_bd(:,2);
    idx = rho_j_r<0;
    rho_j_r(idx) = zeros(size(rho_j_r(idx)))+1;
    
    rho_j_l = (j_bd-0)./e_bd(:,2);
    idx = rho_j_l<0;
    rho_j_l(idx) = zeros(size(rho_j_l(idx)))+1;
    
    
    

    rho_l_bd = min([rho_i_r';rho_i_l';rho_j_r';rho_j_l';rho_l_bd'])';
    
    rho_l(idx_bd) = rho_l_bd;
    
    x = (i_bd - rho_l_bd.*e_bd(:,1))*h;
    y = (j_bd - rho_l_bd.*e_bd(:,2))*h;
    u_l(idx_bd) = g_D(x,y);
    
    
    
    
 %% right   
    % find the inner point
    idx = and(and(1<=i_r,i_r<=N-1),and(1<=j_r,j_r<=N-1));
    s = LocaltoGlobalidx(i_r(idx),j_r(idx),N);
    u_r(idx) = u(s);
    
    idx_bd = ~idx;
    e_bd = e(idx_bd,:);
    
    i_bd = i(idx_bd);
    j_bd = j(idx_bd);
    
     rho_r_bd = rho_r(idx_bd);
    
    rho_i_r = (N-i_bd)./e_bd(:,1);
    idx = rho_i_r<0;
    rho_i_r(idx) = zeros(size(rho_i_r(idx)))+1;
    
    rho_i_l = (0-i_bd)./e_bd(:,1);
    idx = rho_i_l<0;
    rho_i_l(idx) = zeros(size(rho_i_l(idx)))+1;
    
    
   
    rho_j_r = (N-j_bd)./e_bd(:,2);
    idx = rho_j_r<0;
    rho_j_r(idx) = zeros(size(rho_j_r(idx)))+1;
    
    rho_j_l = (0-j_bd)./e_bd(:,2);
    idx = rho_j_l<0;
    rho_j_l(idx) = zeros(size(rho_j_l(idx)))+1;
    
    
    

    rho_r_bd = min([rho_i_r';rho_i_l';rho_j_r';rho_j_l';rho_r_bd'])';
    
    rho_r(idx_bd) = rho_r_bd;

    x = (i_bd+rho_r_bd.*e_bd(:,1))*h;
    y = (j_bd+rho_r_bd.*e_bd(:,2))*h;
    u_r(idx_bd) = g_D(x,y);




Delta_eu = ((u_r-u)./rho_r-(u-u_l)./rho_l).*(2./((rho_l+rho_r).*(e(:,1).^2+e(:,2).^2).*h^2));

if nargout == 2
    Gradu_Deltaeu = [1./rho_l.*(2./((rho_l+rho_r).*(e(:,1).^2+e(:,2).^2).*h^2)),...
                                    1./rho_r.*(2./((rho_l+rho_r).*(e(:,1).^2+e(:,2).^2).*h^2))];
end
end







