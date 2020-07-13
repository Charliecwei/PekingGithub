function u = Newton_method(u,h,G)
global f  alpha Ftype;
N = 1/h;

[i,j] = GlobaltoLocalidx((1:length(u))',N);
% 

switch Ftype
    case 'MAWS'
        [Fu,Grad_Fu] = Grad_MAWS_h_theta(u,h,G);
    case 'MAWS_delta'
        [Fu,Grad_Fu] = Grad_MAWS_h_theta_delta(u,h,G);
end
%cond(full(Grad_Fu))
v = Grad_Fu\(Fu-f(i*h,j*h));



u = u-alpha*v;










%u = u+dt*(MAWS_h_theta(u,h,G)-f(i*h,j*h));
%u = u+dt*(MAWS_h_theta_delta(u,h,G,delta)-f(i*h,j*h));


end