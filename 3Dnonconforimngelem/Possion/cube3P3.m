 function cube3P3
%3-D上P3协调元例子


close all;
format long;
%% Parameters
maxIt = 4; 
N = zeros(maxIt,1); 
h = zeros(maxIt,1);
%% Generate an initial mesh 
[node,elem] = cubemesh([0,1,0,1,0,1],1);

%% Get the data of the pde
%  pde = sincosdata3;
% pde = polydata3;
   pde = Poisson3Data5;
% % %% Set up boundary condition
bdFlag = setboundary3(node,elem,'Dirichlet');




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
    [soln,eqn] = Poisson3P3(node,elem,bdFlag,pde);
    uh = soln.u;
    N(k) = length(uh);
    h(k) = 1./(size(node,1)^(1/3)-1);    
%     errH1(k) = getH1error3(node,elem,pde.Du,uh);
%     errL2(k) = getL2error3(node,elem,pde.exactu,uh);    
%     uI = pde.exactu([node; (node(eqn.edge(:,1),:)+node(eqn.edge(:,2),:))/2]);

    errH1(k) = getH1error3(node,elem,pde.exactDu,uh);
    errL2(k) = getL2error3(node,elem,pde.exactu,uh);    
%     uI =chazhip3erjie(pde.exactu,node,eqn.face,elem);
%     erruIuh(k) = sqrt((uh-uI)'*eqn.Lap*(uh-uI));
      erruIuh(k) =chazhip3wucha(pde.exactu,node,eqn.face,eqn.Lap,eqn.edge,uh);
%     [x,y]=max(abs(uI(1:size(node,1))-uh(1:size(node,1))))
%     size(node,1)
%     size(elem,1)
%     size(eqn.face,1)
    %error(k) = norm(u(1:N(k))-uI(1:N(k)));
     
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
% 
% 
% 
function err=chazhip3wucha(exactu,node,face,Lap,edge,uh)
%计算插值函数误差
%初使参数
%% 计算插值
%点上
un=exactu(node);


%面上插值
uf=exactu((node(face(:,1),:)+node(face(:,2),:)+node(face(:,3),:))/3);

%边上插值
edges=[edge(:,[1,2]);edge(:,[2,1])];

ue=exactu((2*node(edges(:,1),:)+node(edges(:,2),:))/3);

uI=[un;uf;ue];


err=(uI-uh)'*Lap*(uI-uh);
end



function err =getH1error3(node,elem,Du,uh)
%得到P3协调元的一阶误差
%%初使参数 
elem3dof = dof3P3(elem);
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


    la1=lambda(p,1);
    la2=lambda(p,2);
    la3=lambda(p,3);
    la4=lambda(p,4);

    Dla1=Dlambda(:,:,1);
    Dla2=Dlambda(:,:,2);
    Dla3=Dlambda(:,:,3);
    Dla4=Dlambda(:,:,4);

    Dphip = getDphip(Dla1,Dla2,Dla3,Dla4,la1,la2,la3,la4);



       Duh=0;
       for i=1:20
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
elem3dof = dof3P3(elem);

%3阶数值积分
quadOrder = 5;

%% compute L2 error element-wise using quadrature rule with order quadOrder
err = zeros(NT,1);
[lambda,weight] = quadpts3(quadOrder);




phi = getphi(lambda);






 
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
volume = abs(dot(mycross(d12,d13,2),d14,2)/6);
err(isnan(err)) = 0; % singular point is excluded
err = sqrt(sum(volume.*err));
     
                 



end

              
              
function phi = getphi(lambda)

     la1 = lambda(:,1);
     la2 = lambda(:,2);
     la3 = lambda(:,3);
     la4 = lambda(:,4);


% % % phi(:,1) = 9*la1.*(la1-1/3).*(la1-2/3)/2;
% % % phi(:,2) = 9*la2.*(la2-1/3).*(la2-2/3)/2;
% % % phi(:,3) = 9*la3.*(la3-1/3).*(la3-2/3)/2;
% % % phi(:,4) = 9*la4.*(la4-1/3).*(la4-2/3)/2;
% % % 
% % % phi(:,5) = 27*la2.*la3.*la4;
% % % phi(:,6) = 27*la1.*la3.*la4;
% % % phi(:,7) = 27*la1.*la2.*la4;
% % % phi(:,8) = 27*la1.*la2.*la3;
% % % 
% % % phi(:,9) = 27*la1.*la2.*(la1-1/3)/2;
% % % phi(:,10) = 27*la1.*la2.*(la2-1/3)/2;
% % % 
% % % phi(:,11) = 27*la1.*la3.*(la1-1/3)/2;
% % % phi(:,12) = 27*la1.*la3.*(la3-1/3)/2;
% % % 
% % % phi(:,13) = 27*la1.*la4.*(la1-1/3)/2;
% % % phi(:,14) = 27*la1.*la4.*(la4-1/3)/2;
% % % 
% % % phi(:,15) = 27*la2.*la3.*(la2-1/3)/2;
% % % phi(:,16) = 27*la2.*la3.*(la3-1/3)/2;
% % % 
% % % phi(:,17) = 27*la2.*la4.*(la2-1/3)/2;
% % % phi(:,18) = 27*la2.*la4.*(la4-1/3)/2;
% % % 
% % % phi(:,19) = 27*la3.*la4.*(la3-1/3)/2;
% % % phi(:,20) = 27*la3.*la4.*(la4-1/3)/2;



phi(:,1)=dian(la1);
phi(:,2)=dian(la2);
phi(:,3)=dian(la3);
phi(:,4)=dian(la4);

phi(:,5)=mian(la2,la3,la4);
phi(:,6)=mian(la3,la4,la1);
phi(:,7)=mian(la4,la1,la2);
phi(:,8)=mian(la1,la2,la3);

phi(:,[9,10])=bian(la1,la2);
phi(:,[11,12])=bian(la1,la3);
phi(:,[13,14])=bian(la1,la4);
phi(:,[15,16])=bian(la2,la3);
phi(:,[17,18])=bian(la2,la4);
phi(:,[19,20])=bian(la3,la4);



   %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%计算顶点处值
function phis=dian(la)
phis=0.5*la.*(3*la-1).*(3*la-2);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%计算'1，2，3'为面上的基函数值
    function phis=mian(la1,la2,la3)
        phis=27*la1.*la2.*la3;
    end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%计算'1，2'边上基函数值
    function phis=bian(la1,la2)
        phis(:,1)=4.5*la1.*la2.*(3*la1-1);
        phis(:,2)=4.5*la2.*la1.*(3*la2-1);
    end

end
function Dphip = getDphip(Dla1,Dla2,Dla3,Dla4,la1,la2,la3,la4)
% % % 
% % % Dphip(:,:,1) = 9/2*(Da1.*(la1-1/3).*(la1-2/3)+Da1.*la1.*(la1-2/3)+Da1.*la1.*(la1-1/3));
% % % Dphip(:,:,2) = 9/2*(Da2.*(la2-1/3).*(la2-2/3)+Da2.*la2.*(la2-2/3)+Da2.*la2.*(la2-1/3));
% % % Dphip(:,:,3) = 9/2*(Da3.*(la3-1/3).*(la3-2/3)+Da3.*la3.*(la3-2/3)+Da3.*la3.*(la3-1/3));
% % % Dphip(:,:,4) = 9/2*(Da4.*(la4-1/3).*(la4-2/3)+Da4.*la4.*(la4-2/3)+Da4.*la4.*(la4-1/3));
% % % 
% % % Dphip(:,:,5) = 27*(Da2.*la3.*la4+la2.*Da3.*la4+la2.*Da4.*la3);
% % % Dphip(:,:,6) = 27*(Da3.*la4.*la1+la3.*Da4.*la1+la3.*Da4.*la1);
% % % Dphip(:,:,7) = 27*(Da4.*la1.*la2+la4.*Da1.*la2+la4.*Da1.*la2);
% % % Dphip(:,:,8) = 27*(Da1.*la2.*la3+la1.*Da2.*la3+la1.*Da2.*la3);
% % % 
% % % Dphip(:,:,9) = 27/2*(Da1.*la2.*(la1-1/3)+la1.*Da2.*(la1-1/3)+la1.*la2.*Da1);                       
% % % Dphip(:,:,10) = 27/2*(Da1.*la2.*(la2-1/3)+la1.*Da2.*(la2-1/3)+la1.*la2.*Da2);
% % % 
% % % Dphip(:,:,11) = 27/2*(Da1.*la3.*(la1-1/3)+la1.*Da3.*(la1-1/3)+la1.*la3.*Da1);                       
% % % Dphip(:,:,12) = 27/2*(Da1.*la3.*(la3-1/3)+la1.*Da3.*(la3-1/3)+la1.*la3.*Da3);
% % % 
% % % Dphip(:,:,13) = 27/2*(Da1.*la4.*(la1-1/3)+la1.*Da4.*(la1-1/3)+la1.*la4.*Da1);                       
% % % Dphip(:,:,14) = 27/2*(Da1.*la4.*(la4-1/3)+la1.*Da4.*(la4-1/3)+la1.*la4.*Da4);
% % % 
% % % Dphip(:,:,15) = 27/2*(Da2.*la3.*(la2-1/3)+la2.*Da3.*(la2-1/3)+la2.*la3.*Da2);                       
% % % Dphip(:,:,16) = 27/2*(Da2.*la3.*(la3-1/3)+la2.*Da3.*(la3-1/3)+la2.*la3.*Da3);
% % % 
% % % Dphip(:,:,17) = 27/2*(Da2.*la4.*(la2-1/3)+la2.*Da4.*(la2-1/3)+la2.*la4.*Da2);                       
% % % Dphip(:,:,18) = 27/2*(Da2.*la4.*(la4-1/3)+la2.*Da4.*(la4-1/3)+la2.*la4.*Da4);
% % % 
% % % Dphip(:,:,19) = 27/2*(Da3.*la4.*(la3-1/3)+la3.*Da4.*(la3-1/3)+la3.*la4.*Da3);                       
% % % Dphip(:,:,20) = 27/2*(Da3.*la4.*(la4-1/3)+la3.*Da4.*(la4-1/3)+la3.*la4.*Da4);
% % % 
% % % 




Dphip(:,:,1)=dingdian(la1,Dla1);
Dphip(:,:,2)=dingdian(la2,Dla2);
Dphip(:,:,3)=dingdian(la3,Dla3);
Dphip(:,:,4)=dingdian(la4,Dla4);

Dphip(:,:,5)=mian(la2,la3,la4,Dla2,Dla3,Dla4);
Dphip(:,:,6)=mian(la3,la4,la1,Dla3,Dla4,Dla1);
Dphip(:,:,7)=mian(la4,la1,la2,Dla4,Dla1,Dla2);
Dphip(:,:,8)=mian(la1,la2,la3,Dla1,Dla2,Dla3);

Dphip(:,:,[9,10])=bian(la1,la2,Dla1,Dla2);

Dphip(:,:,[11,12])=bian(la1,la3,Dla1,Dla3);

Dphip(:,:,[13,14])=bian(la1,la4,Dla1,Dla4);

Dphip(:,:,[15,16])=bian(la2,la3,Dla2,Dla3);

Dphip(:,:,[17,18])=bian(la2,la4,Dla2,Dla4);

Dphip(:,:,[19,20])=bian(la3,la4,Dla3,Dla4);








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%以la=lambda(:,i)，Dla=Dlambda(:,:,i)算顶点i的基函数导数
function  Dphips=dingdian(la,Dla)
Dphips=0.5*Dla.*(27*power(la,2)-18*la+2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%算'1,2,3'点对应的面上的基函数导数
function Dphips=mian(la1,la2,la3,Dla1,Dla2,Dla3)
   Dphips=27*(Dla1.*la2.*la3+Dla2.*la3.*la1+Dla3.*la1.*la2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%算'1,2'边上的基函数导数
    function Dphips=bian(la1,la2,Dla1,Dla2)
        Dphips(:,:,1)=0.5*Dla2.*(27*power(la1,2)-9*la1)+...
                         0.5*la2.*(54*la1-9).*Dla1;
        Dphips(:,:,2)=0.5*Dla1.*(27*power(la2,2)-9*la2)+...
                         0.5*la1.*(54*la2-9).*Dla2;
    end
end
         
              


