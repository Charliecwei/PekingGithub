function S = Grad_uOfDelta_eu(u,e,h)
%Calculate the gradient _u of Deltaeu (it is a sparse matrix)
[~, Gradu_Deltaeu] = Delta_e(u,e,h);

Nu = length(u);
N = 1/h;

Diag = -sum(Gradu_Deltaeu,2);
Left = Gradu_Deltaeu(:,1);
Right = Gradu_Deltaeu(:,2);
[i,j] = GlobaltoLocalidx((1:Nu)',N);


ii = (1:Nu)';
jj = (1:Nu)';
S = sparse(ii,jj,Diag,Nu,Nu);
%% the left item

    i_l = i- e(:,1);%the x value of the interpolation point on the left
    j_l = j- e(:,2);%the y value of the interpolation point on the left
    
     % find the inner point for the left point
    idx = and(and(1<=i_l,i_l<=N-1),and(1<=j_l,j_l<=N-1));
    jj_left = LocaltoGlobalidx(i_l(idx),j_l(idx),N);
    ii_left = find(idx);
    Left = Left(idx);
    
    S = S + sparse(ii_left,jj_left,Left,Nu,Nu);
%% the right item
    i_r = i+ e(:,1);%the x value of the interpolation point on the left
    j_r = j+ e(:,2);%the y value of the interpolation point on the left
    
    % find the inner point for the right point
    idx = and(and(1<=i_r,i_r<=N-1),and(1<=j_r,j_r<=N-1));
    jj_right = LocaltoGlobalidx(i_r(idx),j_r(idx),N);
    ii_right = find(idx);
    Right = Right(idx);
    
    S = S + sparse(ii_right,jj_right,Right,Nu,Nu);
  %  S = sparse([ii;ii_left;ii_right],[jj;jj_left;jj_right],[Diag;Left;Right],Nu,Nu);
    
end
    
    
    
    





