function example = example1
%% This is the pde data with propagation of a wave packet p69

T = 0.2;
I = [0,1];
h = 0.001;
gama =1.4;
nv = 0.8;


mesh = struct('h',h,'T',T,'I',I,'nv',nv);
solver = {'LF','LW'};
option = struct('solver',{solver});
%option = {'LW'};
pde = struct('uo',@uo,'f',@f,'a',@a,'gama',gama,'bdtype','Inf');

example = struct('mesh',mesh,'option',option,'pde',pde);
  



        
end

function U = uo(x)
%% initial uo
U = zeros(3,size(x,2));
idx = x<0.3;
U(:,idx) = U(:,idx)+[1;0;2.5];
U(:,~idx) = U(:,~idx)+[0.125;0;0.25];
end


function F = f(U,gama)
rhou2 = U(2,:).^2./U(1,:); %rho*u^2;
p = (gama-1)*(U(3,:)-0.5*rhou2);
F = [U(2,:);rhou2+p;(U(3,:)+p).*U(2,:)./U(1,:)];
end


function A = a(U,gama)

A = zeros(3,size(U,2),3); %a(i,j)=A(i,:,j) -- partial(f_i)/partial(u_j);

rho = U(1,:);
u = U(2,:)./rho;
E = U(3,:);

A(2,:,1) = (gama-3)*u.^2/2;
A(3,:,1) = -gama*E.*u./rho+(gama-1)*u.^3;

A(1,:,2) = A(1,:,2)+1;
A(2,:,2) = (3-gama)*u;
A(3,:,2) = gama*E./rho-1.5*(gama-1).*u.^2;

A(2,:,3) = A(2,:,3) + gama-1;
A(3,:,3) = gama*u;



end

