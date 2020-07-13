function DuichenStokes3DP3jia5B4
%%P2jia3B3-P1非协调元求解sotkes方程
% 
%       -div[2*mu*eplison(u)] + grad p = f in \Omega,
%                        - div u = 0  in \Omega,
%   with 
%       Dirichlet boundary condition        u = g_D  on \Gamma_D, 
%       Neumann boundary condition sigma  = g_N  on \Gamma_N. 
% 
%     eplison(u) = 0.5*[grad(u)+grad(u)^T]
%     sigma = (2*mu*eplison(u)- p*I)*n
%%清除数据&加载子文件
close all; 
clear variables;
addpath('../Solver');
addpath('../dof');
addpath('../DuichenStokesData');
addpath('../err');
addpath('../Phi')
%% Parameters
maxIt = 3; 
N = zeros(maxIt,1); 
h = zeros(maxIt,1);
errH1u = zeros(maxIt,1); 
errL2u = zeros(maxIt,1);
erruIh = zeros(maxIt,1);
errp = zeros(maxIt,1);
rateuL2 = zeros(maxIt,1);
rateuH1 = zeros(maxIt,1);
ratepL2 = zeros(maxIt,1);


%% Generate an initial mesh 
%[node,elem] = cubemesh([0,1,0,1,0,1],1);
[node,elem] = cubemesh([-1,1,-1,1,-1,1],2);


%% PDE and options
pde = Duichenstokes3data8;
%pde = stokes3data5;
%% Finite Element Method
for k = 1:maxIt
    % refine mesh
    [node,elem] = uniformrefine3(node,elem);
   
%  bdFlag = setboundary3(node,elem,'Dirichlet');
   bdFlag = setboundary3(node,elem,'Dirichlet','all','Neumann','(0<y&y<1)&(0<x&x<1)&(z==1)');
%   bdFlag = setboundary3(node,elem,'Dirichlet','all','Neumann','z==1'); 
    
    % solve the equation
   [soln,eqn] = DuichenStoke3P3jia5B4andP2(node,elem,bdFlag,pde);

    uh = soln.u;
    ph = soln.p;
    Nu = soln.Nu;
    face = eqn.face;
    N(k) = length(uh)+length(ph);
    h(k) = 2./(power(size(node,1),1/3)-1);
    

    
    
    
     % compute error
     %uI--插值
     uI = chazhi(Nu,pde.exactu,node,face,elem);
    
     %%%误差     
     erruIh(k) = sqrt((uh-uI)'*eqn.A*(uh-uI));
     errp(k) = getL2error3D(node,elem,pde.exactp,ph,'P2');
     errL2u(k) = getstokeserrL2u3D(node,elem,pde.exactu,uh,'P3jia5B4');
     errH1u(k) = getstkoeserrH1u3D(node,elem,pde.exactDu,uh,'P3jia5B4');
    
    
end


for k = 2:maxIt
    rateuL2(k) = (log(errL2u(k-1))-log(errL2u(k)))/log(2);
    ratepL2(k) = (log(errp(k-1))-log(errp(k)))/log(2);
    rateuH1(k) = (log(errH1u(k-1))-log(errH1u(k)))/log(2);
end



huitu(h,errp,errL2u,errH1u,rateuL2,rateuH1,ratepL2, N);

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%辅助函数

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%绘图
function huitu(h,errp,errL2u,errH1u,rateuL2,rateuH1,ratepL2, N)
%%%% Plot convergence rates


showrateh3(h,errH1u,2,'-*', '| u - u_h |_1', ...
           h,errL2u,2,'k-+', '|| u - u_h ||', ...
           h,errp,2,'m-+','|| p - p_h ||');

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','|| u-u_h ||','rateuL2','| u-u_h |_1','rateuH1','|| p-p_h ||','ratepL2'};
disptable(colname,N,[],h,'%0.3e',errL2u,'%0.3e',rateuL2,'%0.2f',errH1u,'%0.3e',rateuH1,'%0.2f',...
              errp,'%0.3e',ratepL2,'%0.2f');     
%saveas(gcf,'dane6','png')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%插值uI
function uI = chazhi(Nu,exactu,node,face,elem)
      NQ = Nu-size(elem,1);
      u1 = zeros(Nu,1);
      u2 = zeros(Nu,1);
      u3 = zeros(Nu,1);
     [lambda,w] = quadpts(5);
     nQuad=length(lambda);
     
     for p=1:nQuad
            pxyz=lambda(p,1)*node(face(:,1),:)+lambda(p,2)*node(face(:,2),:)...
                                +lambda(p,3)*node(face(:,3),:);
             uD=exactu(pxyz);
             uD1=uD(:,1);
             uD2=uD(:,2);
             uD3=uD(:,3);
             
             la1 = lambda(p,1);
             la2 = lambda(p,2);
             la3 = lambda(p,3);
             
             u1(1:NQ)=u1(1:NQ)+w(p)*[uD1*la1^2;uD1*la2^2;uD1*la3^2;uD1*la1*la2;uD1*la1*la3;uD1*la2*la3];
             u2(1:NQ)=u2(1:NQ)+w(p)*[uD2*la1^2;uD2*la2^2;uD2*la3^2;uD2*la1*la2;uD2*la1*la3;uD2*la2*la3];
             u3(1:NQ)=u3(1:NQ)+w(p)*[uD3*la1^2;uD3*la2^2;uD3*la3^2;uD3*la1*la2;uD3*la1*la3;uD3*la2*la3];
 
     end
     
         [lambda,w] = quadpts3(5);
         nQuad=length(lambda);

       for p=1:nQuad
          pxyz=lambda(p,1)*node(elem(:,1),:)+lambda(p,2)*node(elem(:,2),:)...
                               +lambda(p,3)*node(elem(:,3),:)+lambda(p,4)*node(elem(:,4),:);
             uD=exactu(pxyz);
             uD1=uD(:,1);
             uD2=uD(:,2);
             uD3=uD(:,3);
             
             u1(NQ+1:end)=u1(NQ+1:end)+w(p)*uD1;
             u2(NQ+1:end)=u2(NQ+1:end)+w(p)*uD2;
             u3(NQ+1:end)=u3(NQ+1:end)+w(p)*uD3;
             
      end
      uI = [u1;u2;u3];

end


