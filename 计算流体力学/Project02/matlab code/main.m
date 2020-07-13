%example = example1;
example = example2;








i = 2;

pde = example.pde;
option = example.option;
mesh = example.mesh;

options = struct('solver',option.solver{i});



[x,U] = ConservativeFormat(pde,options,mesh);
rho = U(1,:);
u = U(2,:)./rho;
E = U(3,:);
p = (pde.gama-1)*(E-0.5*rho.*u.^2);

plot(x,rho);hold on
plot(x,u);hold on
plot(x,p);hold on

%plot(mesh.I(1):mesh.h:mesh.I(2),E,'--');hold on
legend(['rho, CFL = ',num2str(mesh.nv)],['u, CFL = ',num2str(mesh.nv)],['p, CFL = ',num2str(mesh.nv)]);
title([option.solver{i},', h is ',num2str(mesh.h)])

% 
%  legend('CFL = 0.6','CFL = 0.8','CFL = 1');
% title(['u, ',option.solver{i},', h is ',num2str(mesh.h)])





