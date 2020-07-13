function cube3DP3N
%2-D上P2非协调元例子

addpath('../dof');
addpath('../Phi');
addpath('../PossionData');
addpath('../Solver');
addpath('../err');
close all;
format long;
%% Parameters
maxIt = 3; 
N = zeros(maxIt,1); 
h = zeros(maxIt,1);

%% Generate an initial mesh 
[node,elem] = cubemesh([0,1,0,1,0,1],1);

%% Get the data of he pde
pde = Poisson3Data1;
% pde = polydata3;


%bdFlag = setboundary3(node,elem,'Dirichlet');
%bdFlag = setboundary3(node,elem,'Neumann');
bdFlag = setboundary3(node,elem,'Dirichlet','all','Neumann','(x==1|y==1|z==1)');

   

%% Finite Element Method        
errL2 = zeros(maxIt,1); 
errH1 = zeros(maxIt,1); 
rateuL2 = zeros(maxIt,1);
rateuH1 = zeros(maxIt,1);




for k=1:maxIt
    % refine grid    
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag); 
    %% Set up boundary condition

    % solve the equation
    [soln,eqn] = Poisson3DP3N(node,elem,bdFlag,pde); 
    uh = soln.u;
    N(k) = length(uh);
    h(k) = 1./(size(node,1)^(1/3)-1);    

    errH1(k) = getH1error3D(node,elem,pde.exactDu,uh,'3DP3N');
    errL2(k) = getL2error3D(node,elem,pde.exactu,uh,'3DP3N');    
  
 
end

%% Plot convergence rates

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
