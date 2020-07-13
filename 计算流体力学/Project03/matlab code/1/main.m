example = example1;
%example = example2;



global solver phivar;
solver = 'LF';% solver = 'LF';
%phi = { 'zero','one','minmod','Superbee'};
% phi = { 'zero', 'minmod','Superbee'};
phi = { 'zero','one','minmod','Superbee'};








pde = example.pde;
mesh = example.mesh;





%%

S = cell(1,length(phi));
for i = 1:length(phi)
    phivar = phi{i};
    [x,U] = ConservativeFormat(pde,mesh);
    rho1 = U(1,:);
    u1 = U(2,:)./rho1;
    E1 = U(3,:);
    p1 = (pde.gama-1)*(E1-0.5*rho1.*u1.^2);
    figure(1)
    plot(x,rho1,'--');hold on
    figure(2)
    plot(x,u1,'--');hold on
    figure(3)
    plot(x,p1,'--');hold on
    S{i} = ['phi is ',phi{i}];
end


% %% exact solution
% U = pde.uo([0.2,0.4]);
% rho = U(1,:);
% u = U(2,:)./rho;
% E = U(3,:);
% p = (pde.gama-1)*(E-0.5*rho.*u.^2);
% 
% tol = 1e-5;
% rho1 = rho(1);
% rho4 = rho(end);
% u1 = u(1);
% u4 = u(end);
% p1 = p(1);
% p4 = p(end);
% 
% x0 = -0.3;
% xf = 0.7;
% t = mesh.T;
% npoints = 1/mesh.h;
% [ERho,EU,EP] = RiemannExact(p1,rho1,u1,p4,rho4,u4,tol,t,x0,xf,npoints);
% figure(1)
% plot(x,ERho,'-b','LineWidth',2);
% figure(2)
% plot(x,EU,'-b','LineWidth',2);
% figure(3)
% plot(x,EP,'-b','LineWidth',2);
% S{length(phi)+1} = 'exact';
% 



figure(1)
legend(S);

figure(2)
legend(S);

figure(3)
legend(S);



figure(1)
title(['rho, h is ',num2str(mesh.h)]);

figure(2)
title(['u, h is ',num2str(mesh.h)]);

figure(3)
title(['p, h is ',num2str(mesh.h)]);









