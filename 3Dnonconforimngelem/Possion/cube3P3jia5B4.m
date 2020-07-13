 function cube3P3jia5B4
%3-D上P3协调元例子
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

%% Get the data of the pde
%  pde = sincosdata3;
% pde = polydata3;
  % pde = x5eycosz3;
  pde = Poisson3Data5;
% % %% Set up boundary condition
bdFlag = setboundary3(node,elem,'Dirichlet');



%% Finite Element Method  

errL2 = zeros(maxIt,1); 
errH1 = zeros(maxIt,1); 
erruIuh = zeros(maxIt,1);

for k=1:maxIt
    % refine grid    
    
    % solve the equation
    [soln,eqn] = Poisson3P3jia5B4(node,elem,bdFlag,pde);
    uh = soln.u;
    N(k) = length(uh);
    h(k) = 1./(size(node,1)^(1/3)-1);    

    errH1(k) = getH1error3(node,elem,pde.exactDu,uh);
    errL2(k) = getL2error3(node,elem,pde.exactu,uh);    

    erruIuh(k) =chazhip3wucha(pde.exactu,node,eqn.face,eqn.Lap,uh,elem);
    [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);

     
 end

%%%%%%%%Plot convergence rates
figure;
showrateh2(h,errH1,1,'-*', '|| Du - Du_h ||', ...
           h,errL2,1,'k-+', '|| u - u_h ||');%, ...
          % h,erruIuh,1,'m-+','|| D u_I - D u_h ||');

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','|| u-u_h ||','|| Du-Du_h ||'};%,'|| Du_I-Du_h ||'};
disptable(colname,N,[],h,'%0.3e',errL2,'%0.5e',errH1,'%0.5e')%,erruIuh,'%0.5e');

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
       ge = ge.*repmat(volume,1,1);
       uI(NQ+1:end) = ge(:);

err=(uI-uh)'*Lap*(uI-uh);
end

function err =getH1error3(node,elem,Du,uh)
%得到P3协调元的一阶误差
%%初使参数 
elem3dof = dof3P3jia5B4s(elem);
NT = size(elem,1); 
%3阶数值积分
quadOrder = 5;
%计算梯度和体积
[Dlambda,volume] = gradbasis3(node,elem);


%% compute H1 error element-wise using quadrature rule with order quadOrder
[lambda,weight] = quadpts3(quadOrder);
nQuad = size(lambda,1);
err = zeros(NT,1);
for p = 1:nQuad
    % quadrature points in the x-y-z coordinate
    pxyz = lambda(p,1)*node(elem(:,1),:) ...
         + lambda(p,2)*node(elem(:,2),:) ...
         + lambda(p,3)*node(elem(:,3),:) ...
         + lambda(p,4)*node(elem(:,4),:);


%     la1=lambda(p,1);
%     la2=lambda(p,2);
%     la3=lambda(p,3);
%     la4=lambda(p,4);
% 
%     Dla1=Dlambda(:,:,1);
%     Dla2=Dlambda(:,:,2);
%     Dla3=Dlambda(:,:,3);
%     Dla4=Dlambda(:,:,4);

    Dphip = getDphipP3jiapos(Dlambda,lambda(p,:));



       Duh=0;
       for i=1:25
             Duh =Duh+ repmat(uh(elem3dof(:,i)),1,3).*Dphip(:,:,i);
       end
              
              
       
       
       err = err + weight(p)*sum(((Du(pxyz)-Duh).^2),2);
end

err = volume.*err;
err(isnan(err)) = 0; % singular values are excluded

err = sqrt(sum(err));




end

function err=getL2error3(node,elem,exactu,uh)
%得到P3三次协调元的零阶误差
%%初使参数 
NT = size(elem,1); 
elem3dof = dof3P3jia5B4s(elem);

%3阶数值积分
quadOrder = 5;

%% compute L2 error element-wise using quadrature rule with order quadOrder
err = zeros(NT,1);
[lambda,weight] = quadpts3(quadOrder);




phi = getphiP3jiapos(lambda);






 
 nQuad = size(lambda,1);
for p = 1:nQuad
    % evaluate uh at quadrature point
     uhp = sum(uh(elem3dof(:,:)).*phi(p,:) ,2);
     % quadrature points in the x-y coordinate
    pxy = lambda(p,1)*node(elem(:,1),:) ...
        + lambda(p,2)*node(elem(:,2),:) ...
        + lambda(p,3)*node(elem(:,3),:) ...
        + lambda(p,4)*node(elem(:,4),:);
    err = err + weight(p)*((exactu(pxy) - uhp).^2);
end
%% Modification
% volume of tetrahedrons
d12 = node(elem(:,2),:)-node(elem(:,1),:);
d13 = node(elem(:,3),:)-node(elem(:,1),:);
d14 = node(elem(:,4),:)-node(elem(:,1),:);
volume = abs(dot(mycross(d12,d13,2),d14,2))/6;
err(isnan(err)) = 0; % singular point is excluded
err = sqrt(sum(volume.*err));
     
                 



end             
