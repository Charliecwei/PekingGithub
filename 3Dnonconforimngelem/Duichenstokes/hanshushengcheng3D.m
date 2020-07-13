clear
syms x y z real;
%u=curl(u);
X = [x,y,z];
%u = [ x.^3.*y*pi.*sin(pi*z),3*z.*(x - 1).^2 + y.^3.*exp(x),3*x.^2.*y.*cos(pi*z) - 3*y.^2.*z.*exp(x)];
u = [y.^6.*(1-y).^6*x.^5.*(1-x).^5.*z^3.*(1-z).^4,y.^5.*(1-y).^5*x.^6.*(1-x).^6.*z^3.*(1-z).^4,0];
u = curl(u)';
p=(x-1/2)*(y-1/2)*(1-z);


Du = grad(u,X);
Eu = (Du+Du')/2;
DivDu = divce(Du,X);
DivEu = divce(Eu,X);
Dp = grad(p,X)';

f = -DivDu+Dp';
Duichengf = -2*DivEu+Dp';

divce(u,X);
%%
%求u的梯度. Du(i,j) = parital_i(u_j)
function Du = grad(u,X)
n = length(X);
m = size(u,2);
%Du = zeros(n);
%syms Du real;
for i = 1:n
    for j = 1:m
        Du(i,j) = diff(u(j),X(i));
    end
end
Du = Du';
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
Divu= Divu';
end

