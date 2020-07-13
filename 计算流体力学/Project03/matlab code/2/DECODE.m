function [u,v,p,T] = DECODE(U1,U2,U3,U5)
%% git rho u v e from U1 U2 U3 U4
global c_v R;

rho = U1;
u = U2./rho;
v = U3./rho;

e = U5./rho - 0.5*(u.^2+v.^2);

T = e./c_v;
p = rho.*R.*T;
% sum(e(:)<0)

end