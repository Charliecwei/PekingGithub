clear
syms x y z real;
%u=curl(u);
X = [x,y,z];
% u = [ x.^3.*y*pi.*sin(pi*z),3*z.*(x - 1).^2 + y.^3.*exp(x),3*x.^2.*y.*cos(pi*z) - 3*y.^2.*z.*exp(x)];
% u = curl(u)';
% p=(x-1/2)*(y-1/2)*(1-z);

% u = [pi*cos(pi*y)*(sin(pi*x)^2)*sin(pi*y)*sin(pi*z),-pi*cos(pi*x)*(sin(pi*y)^2)*sin(pi*x)*sin(pi*z),0];
% p = 0*(x+y+z);

% u = [10*x*y^4+10*x*z^4-4*x^5,10*y*x^4+10*y*z^4-4*y^5,10*z*x^4+10*z*y^4-4*z^5];
% p = -60*x^2*y^2-60*x^2*z^2-60*y^2*z^2+20*(x^4+y^4+z^4);

% u = [2*y*z*(x^2-1)^2*(y^2-1)*(z^2-1),-x*z*(x^2-1)*(y^2-1)^2*(z^2-1),-x*y*(x^2-1)*(y^2-1)*(z^2-1)^2];
% p = x*y*z;

u = [0*x,0*y,0*z];
p = -60*x^2*y^2-60*x^2*z^2-60*y^2*z^2+20*(x^4+y^4+z^4)+8;


Du = grad(u,X);
Eu = (Du+Du')/2;
DivDu = divce(Du,X);
DivEu = divce(Eu,X);
Dp = grad(p,X);

f = -DivDu+Dp;
Duichengf = -2*DivEu+Dp;

divce(u,X)
%%
%求u的梯度. Du(i,j) = parital_i(u_j)
function Du = grad(u,X)
n = length(X);
m = size(u,2);
%Du = zeros(n);
%syms Du real;
for i = 1:m
    for j = 1:n
        Du(i,j) = diff(u(i),X(j));
    end
end
end

%求u的散度.u(i,j)表示第j个函数的第i个分量
function Divu = divce(u,X)
n = length(X);
m = size(u,1);
%syms Divu real;
%Divu = 0;
for i = 1:m
    for j = 1:n
        Divu(i,j) = diff(u(i,j),X(j));
    end
end
Divu = sum(Divu,2);
Divu = Divu';
end

