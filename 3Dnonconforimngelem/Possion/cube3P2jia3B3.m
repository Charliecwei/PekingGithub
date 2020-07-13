%function Poisson3DP2jia
%%%%%%%% 解Poisson方程 3D P2jia元
%%清除数据
close all; 
clear variables;
addpath('../Solver');
addpath('../dof');
addpath('../PossionData');
addpath('../err');
addpath('../Phi');
%% Parameters
maxIt = 4; 
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
[node,elem] = cubemesh([0,1,0,1,0,1],1);


%% PDE and options
pde = Poisson3Data5;
%% Finite Element Method
for k = 1:maxIt
    % refine mesh
   [node,elem] = uniformrefine3(node,elem);
    %[node,elem] = cubemesh([0,1,0,1,0,1],0.5);
    
 bdFlag = setboundary3(node,elem,'Dirichlet');
% bdFlag = setboundary3(node,elem,'Neumann');%--------其相容性条件不好提！！！
%bdFlag = setboundary3(node,elem,'Neumann','all','Dirichlet','y==0');
   
    
    % solve the equation
  [soln,eqn] = Poisson3P2jia3B3(node,elem,bdFlag,pde);

    uh = soln.u;
    Nu = soln.Nu;
    face = eqn.face;
    N(k) = length(uh);
    h(k) = 1./(power(size(node,1),1/3)-1);
    

    
    
    
     % compute error
     %uI--插值
     uI = chazhi(Nu,pde.exactu,node,face,elem);
    
     %%%误差     
     erruIh(k) = sqrt((uh-uI)'*eqn.A*(uh-uI));
     errL2u(k) = getL2error3D(node,elem,pde.exactu,uh,'P2jia');
     errH1u(k) = getH1error3D(node,elem,pde.exactDu,uh,'P2jia');
    
     
end


for k = 2:maxIt
    rateuL2(k) = (log(errL2u(k-1))-log(errL2u(k)))/log(2);
    ratepL2(k) = (log(errp(k-1))-log(errp(k)))/log(2);
    rateuH1(k) = (log(errH1u(k-1))-log(errH1u(k)))/log(2);
end

huitu(h,erruIh,errL2u,errH1u,rateuL2,rateuH1,N);




%end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%辅助函数

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%绘图
function huitu(h,erruIh,errL2u,errH1u,rateuL2,rateuH1,N)
%%%% Plot convergence rates

figure(2);
showrateh3(h,errL2u,1,'m-+','||u-u_h||_2',...
                      h,erruIh,1,'k-+','||I_hu-u_h||_2',...
                      h,errH1u,1,'-*','|u-u_h|_1');
                  
                  
disp('Table: Error')
colname = {'#Dof','h','||u-u_h||_2','rateuL2',...
                                  '|u-uh|_1','rateuH1' ,...
                                  };
disptable(colname,N,[],h,'%0.3e',...
                                  errL2u,'%0.3e',rateuL2,'%0.2f',...
                                   errH1u,'%0.3e',rateuH1,'%0.2f' );
%saveas(gcf,'dane6','png')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%插值uI
function uI = chazhi(Nu,exactu,node,face,elem)
      NQ = Nu-size(elem,1);
      u = zeros(Nu,1);
     [lambda,w] = quadpts(5);
     nQuad=length(lambda);
     
     for p=1:nQuad
            pxyz=lambda(p,1)*node(face(:,1),:)+lambda(p,2)*node(face(:,2),:)...
                                +lambda(p,3)*node(face(:,3),:);
             uD=exactu(pxyz);
             
                          
             u(1:NQ)=u(1:NQ)+w(p)*[uD*lambda(p,1);uD*lambda(p,2);uD*lambda(p,3)];
  
     end
     
         [lambda,w] = quadpts3(5);
         nQuad=length(lambda);

       for p=1:nQuad
          pxyz=lambda(p,1)*node(elem(:,1),:)+lambda(p,2)*node(elem(:,2),:)...
                               +lambda(p,3)*node(elem(:,3),:)+lambda(p,4)*node(elem(:,4),:);
             uD=exactu(pxyz);
             
             u(NQ+1:end)=u(NQ+1:end)+w(p)*uD;
             
      end
      uI = u;

end











