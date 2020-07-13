 function cube3P3jia5B4s
%3-D上P3协调元例子
addpath('../Phi');
addpath('../PossionData');
addpath('../Solver');
addpath('../err');
addpath('../dof');
close all;
format long;
%% Parameters
maxIt = 3; 
N = zeros(maxIt,1); 
h = zeros(maxIt,1);
%% Generate an initial mesh 
[node,elem] = cubemesh([0,1,0,1,0,1],1);

%% Get the data of the pde
%  pde = sincosdata3;
% pde = polydata3;
%   pde = x5eycosz3;
 pde = Poisson3Data5;

%%%%%%% 非齐次的会掉两阶，纯Neumann的不会做!!!
% % %% Set up boundary condition
bdFlag = setboundary3(node,elem,'Dirichlet');
%bdFlag = setboundary3(node,elem,'Neumann','all','Dirichlet','z==0|x==0');
   


%% Finite Element Method  

errL2 = zeros(maxIt,1); 
errH1 = zeros(maxIt,1); 
erruIuh = zeros(maxIt,1);
rateuL2 = zeros(maxIt,1);
rateuH1 = zeros(maxIt,1);


for k=1:maxIt
    % refine grid    
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
    % solve the equation
    [soln,eqn] = Poisson3P3jia5B4(node,elem,bdFlag,pde);
    uh = soln.u;
    N(k) = length(uh);
    h(k) = 1./(size(node,1)^(1/3)-1);    

    errH1(k) = getH1error3D(node,elem,pde.exactDu,uh,'P3jia5B4s');
    errL2(k) = getL2error3D(node,elem,pde.exactu,uh,'P3jia5B4s');    

      erruIuh(k) =chazhip3wucha(pde.exactu,node,eqn.face,eqn.Lap,uh,elem);

      
end
 

for k = 2:maxIt
    rateuL2(k) = (log(errL2(k-1))-log(errL2(k)))/log(2);
    rateuH1(k) = (log(errH1(k-1))-log(errH1(k)))/log(2);
end

%%%%%%%%Plot convergence rates
figure;
showrateh2(h,errH1,1,'-*', '|| Du - Du_h ||', ...
           h,errL2,1,'k-+', '|| u - u_h ||');%, ...
          % h,erruIuh,1,'m-+','|| D u_I - D u_h ||');

disp('Table: Error')
colname = {'#Dof','h','||u-u_h||_2','rateuL2',...
                                  '|u-uh|_1','rateuH1' ,...
                                  };
disptable(colname,N,[],h,'%0.3e',...
                                  errL2,'%0.3e',rateuL2,'%0.2f',...
                                   errH1,'%0.3e',rateuH1,'%0.2f' );

end
% 

function err=chazhip3wucha(exactu,node,face,Lap,uh,elem)
%计算插值函数误差
%初使参数
%% 计算插值
%%面上
[lambda,w] = quadpts(5);
nQuad=length(lambda);
NQ = 6*size(face,1);
Nu = size(uh,1);
uI = zeros(Nu,1);

    for pp=1:nQuad
            pxyz=lambda(pp,1)*node(face(:,1),:)+lambda(pp,2)*node(face(:,2),:)...
                                +lambda(pp,3)*node(face(:,3),:);
             uD=exactu(pxyz);          
             la1 = lambda(pp,1);
             la2 = lambda(pp,2);
             la3 = lambda(pp,3);
             uI(1:NQ)=uI(1:NQ)+w(pp)*[uD*la1^2;uD*la2^2;uD*la3^2;uD*la1*la2;uD*la1*la3;uD*la2*la3];
    end
     
     
    
 %%体上   
 [lambda,w] = quadpts3(5);
  nQuad=length(lambda);
         
         
         ge = 0;

       for pp=1:nQuad
          pxyz=lambda(pp,1)*node(elem(:,1),:)+lambda(pp,2)*node(elem(:,2),:)...
                               +lambda(pp,3)*node(elem(:,3),:)+lambda(pp,4)*node(elem(:,4),:);
             uD=exactu(pxyz);                         
             ge = ge + w(pp)*uD;  
       end
       [~,volume] = gradbasis3(node,elem);
       ge = ge.*volume;
       uI(NQ+1:end) = ge;

err=(uI-uh)'*Lap*(uI-uh);
end





     


