function s = Delta_e(v,e)
%% compute the delta_ev for v on boundary=0

global n N0



    s = -2*v;
    I = (1:N0)';
    j =  ceil(I/(n-1));
    i = I - (j-1)*(n-1);

    % plus  
    ip = i+e(1);
    jp = j+e(2);
    
    idx =  logical((ip>0).*(ip<n).*(jp>0).*(jp<n));
    Np = (jp(idx)-1)*(n-1)+ip(idx);
    s(idx) = s(idx)+v(Np);

    % mins  
    ip = i-e(1);
    jp = j-e(2);
    idx =  logical((ip>0).*(ip<n).*(jp>0).*(jp<n));
    Np = (jp(idx)-1)*(n-1)+ip(idx);
    s(idx) = s(idx)+v(Np);
    
    s = s*(n/2)^2/sum(e.^2);
    
end    
    