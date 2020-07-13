%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%2019������Ԫ�ϻ���ҵһ
close all; 
clear variables;
%%var = 'LagrangeP1','HermiteP3','HermiteP4';��������ʲôԪ
var = 'HermiteP4';
Hw(var);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%������
function Hw(var)



%% Parameters 
maxIt = 4; 
N = zeros(maxIt,1);
h = zeros(maxIt,1);
errL2 = zeros(maxIt,1); 
errH1 = zeros(maxIt,1); 
errLinf = zeros(maxIt,1);
rateL2 = zeros(maxIt,1);
rateH1 = zeros(maxIt,1);
rateLinf = zeros(maxIt,1);
computingtime = zeros(maxIt,1);

%% Generate an initial mesh 
switch var
    case 'LagrangeP1'
            [node,elem] = squaremesh([0 1 0 1], 0.25);
            Sk = 2;
    case 'HermiteP3'
        [node,elem] = squaremesh([0 1 0 1], 0.25);
        Sk = 2;
    case 'HermiteP4'
        [node,elem] = squaremesh([0 1 0 1], 0.5);
        Sk = 1;
end

bdFlag = setboundary(node,elem,'Dirichlet','x==0|x==1','Neumann','y==0|y==1');
for k = 1:Sk
    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
end

%% Get the data of the pde
pde = Possion2data2;


%% Finite Element Method        
for k = 1:maxIt
    % refine mesh
    time = cputime;
   [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
    % solve the equation
    
    switch var
        
        case 'LagrangeP1'     %Lagrange P1 Ԫ
                [soln,eqn,info] = PoissonP1(node,elem,bdFlag,pde); 
     
        case 'HermiteP3'            % HermiteP3Ԫ
               [soln,eqn,info,area] = PoissonHermiteP3s(node,elem,bdFlag,pde);
      
        case 'HermiteP4'         % HermiteP4Ԫ
              [soln,eqn,info,area] = PoissonHermiteP4s(node,elem,bdFlag,pde);
              
    end
    computingtime(k) = cputime-time;
 
    uh = soln.u;
    N(k) = max(double(eqn.elemdof(:)));
    h(k) = 1./(sqrt(size(node,1))-1);
    % compute error
    
    switch var
        
        case 'LagrangeP1'   % LagrangeP1Ԫ���
                errLinf(k) = getLinferror(node,elem,pde.exactu,uh); 
                errL2(k) = getL2error(node,elem,pde.exactu,uh);
                errH1(k) = getH1error(node,elem,pde.exactDu,soln.Du);

        case 'HermiteP3'   %HermiterP3Ԫ���
                errLinf(k) = getLinfH3errors(node,eqn.elemdof,elem,pde.exactu,uh);
                errL2(k) = getL2H3errors(node,eqn.elemdof,elem,area,pde.exactu,uh);
                errH1(k) = getH1H3errors(node,eqn.elemdof,elem,pde.exactDu,uh);

        case 'HermiteP4'    %HermiterP4Ԫ���
                errLinf(k) = getLinfH4errors(node,eqn.elemdof,elem,pde.exactu,uh);
                errL2(k) = getL2H4errors(node,eqn.elemdof,elem,area,pde.exactu,uh);
                errH1(k) = getH1H4errors(node,eqn.elemdof,elem,pde.exactDu,uh);
                
    end
end


%% ��������
for k = 2:maxIt
    rateL2(k) = (log(errL2(k-1))-log(errL2(k)))/log(2);
    rateLinf(k) = (log(errLinf(k-1))-log(errLinf(k)))/log(2);
    rateH1(k) = (log(errH1(k-1))-log(errH1(k)))/log(2);
end



figure(2);
showrateh3(h,errL2,1,'m-+','||u-u_h||_{L2}',...
                      h,errLinf,1,'k-+','||u-uh||_{inf}',...
                      h,errH1,1,'-*','|u-uh|_{H1}');

fprintf('\n');
disp('Table: Error')
colname = {'#Dof','h','cuptime','||u-u_h||_{inf}','rateLinf',...
                                          '||u-uh||_{L2}','rateL2',...
                                           '|u-uh|_{H1}','rateH1'};
disptable(colname,N,[],h,'%0.5e',computingtime,'%0.5e',...
                                                                  errLinf,'%0.5e',rateLinf,'%0.5e',...
                                                                   errL2,'%0.5e',rateL2,'%0.5e', ...
                                                                   errH1,'%0.5e',rateH1,'%0.5e');

end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LagrangeP1
function [soln,eqn,info] = PoissonP1(node,elem,bdFlag,pde,option)
%%2019������Ԫ�������ϻ���ҵ1,ΪP1Ԫ����[0,1]^2�ϵĻ�ϱ߽�������Possion����.

N = size(node,1); %�����;
NT = size(elem,1);%��Ԫ����;
Ndof = N;%ȫ�ֻ���������;
elemdof = elem; %�ֲ���������Ӧ�������������λ��.��elemdof(i,j)����ʾ��i����Ԫ�ĵ�j����������Ӧ��������������;

%%
%%��װ�նȾ���.
%1.�����ݶ������
%��ʱ����ת���,��ev1=(a,b),�����Ӧ��ָ��a1�Ĵ�ֱ����Ϊe1^{T}=(-b,a);���������Ϊarea=e1^{T}*e2/2;
%��grad(lambda_1)=e1^{T}/2area;
ev1 = node(elem(:,3),:)-node(elem(:,2),:);
ev2 = node(elem(:,1),:)-node(elem(:,3),:);
ev3 = node(elem(:,2),:)-node(elem(:,1),:);
area = (-ev1(:,2).*ev2(:,1)+ev1(:,1).*ev2(:,2));

Dlambda(1:NT,:,1) = [-ev1(:,2)./area,ev1(:,1)./area];
Dlambda(1:NT,:,2) = [-ev2(:,2)./area,ev2(:,1)./area];
Dlambda(1:NT,:,3) = [-ev3(:,2)./area,ev3(:,1)./area];
Dphi= Dlambda;
area = 0.5*area;
is = (area<0);
area(is) = -area(is);
%Dlambda(i,:,s)��ʾ��i����Ԫ��lambda_s���ݶȵ���

%�ֲ��նȾ����Լ����Ӧ������նȾ�������
time = cputime;
ii = zeros(6*NT,1);
jj = zeros(6*NT,1);
Aij = zeros(6*NT,1);
Np = 0;
for i = 1:3
    for j = i:3
        ii(Np+1:Np+NT) = elemdof(:,i);
        jj(Np+1:Np+NT) = elemdof(:,j);
        Aij(Np+1:Np+NT) = dot(Dlambda(:,:,i),Dlambda(:,:,j),2).*area;
        Np = Np + NT;
    end
end

%��װ����նȾ���
idx = (ii==jj);
AD = sparse(ii(idx),jj(idx),Aij(idx),Ndof,Ndof);
AU = sparse(ii(~idx),jj(~idx),Aij(~idx),Ndof,Ndof);
A = AD+AU+AU';
clear Aij AD AU
%%


%%�����Ҷ���
option.fquadorder = 3; %���������
[lambda,weight] = quadpts(option.fquadorder);
phi = lambda;%�ֲ�������
nQuad = size(lambda,1);
bt = zeros(NT,3);
 
for pp = 1:nQuad
    pxy = lambda(pp,1).*node(elem(:,1),:)+lambda(pp,2).*...
                    node(elem(:,2),:)+lambda(pp,3).*node(elem(:,3),:);
     fp = pde.f(pxy);
     for i = 1:3
         bt(:,i) = bt(:,i) + weight(pp)*phi(pp,i)*fp;
     end
end
bt = bt.*repmat(area,1,3);
b = accumarray(elemdof(:),bt(:),[Ndof 1]);
clear pxy bt

%%
%%b�߽����.
[AD,b,u,freeNode,fixedNode] = getbd(b);


%%
%%�����Է���
assembleTime = cputime - time;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Solve the system of linear equations
if isempty(freeNode), return; end
% Set up solver type
if isempty(option) || ~isfield(option,'solver')  || isfield(option,'mgoption')   % no option.solver
    if Ndof <= 2e3  % Direct solver for small size systems
        option.solver = 'direct';
    else            % MGCG  solver for large size systems
        option.solver = 'mg';
    end
end
solver = option.solver;
% solve
switch solver
    case 'direct'
        t = cputime;
        u(freeNode) = AD(freeNode,freeNode)\b(freeNode);
        residual = norm(b - AD*u);
        info = struct('solverTime',cputime - t,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
    case 'none'
        info = struct('solverTime',[],'itStep',0,'err',[],'flag',3,'stopErr',[]);
    case 'mg'
        if ~isfield(option,'mgoption')   % no option.mgoption
            option.mgoption.x0 = u;
            option.mgoption.solver = 'CG';
        end
        [u,info] = mg(AD,b,elem,option.mgoption);
    case 'amg'
        if ~isfield(option,'amgoption')  % no option.amgoption
            option.amgoption.x0 = u;
            option.amgoption.solver = 'CG';
        end
        [u(freeNode),info] = amg(AD(freeNode,freeNode),b(freeNode),option.amgoption);                 
end


dudx =  u(elem(:,1)).*Dphi(:,1,1) + u(elem(:,2)).*Dphi(:,1,2) ...
      + u(elem(:,3)).*Dphi(:,1,3);
dudy =  u(elem(:,1)).*Dphi(:,2,1) + u(elem(:,2)).*Dphi(:,2,2) ...
      + u(elem(:,3)).*Dphi(:,2,3);         
Du = [dudx, dudy];

%% Output
if nargout == 1
    soln = u;
else
    soln = struct('u',u,'Du',Du);
    eqn = struct('A',AD,'b',b,'freeNode',freeNode,'elemdof',elemdof,...
                              'fixedNode',fixedNode,'Lap',A);
    info.assembleTime = assembleTime;
end








%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [AD,b,u,freeNode,fixedNode] = getbd(b)
    %%����Dirichle�߽�������Neumann�߽�����

    %1. ��Dirichlet�߽�,����Ϊ:fixedNode,����������A
    [fixedNode,bdEdge,isBdNode] = findboundary(elem,bdFlag);
    freeNode = find(~isBdNode);


    %�ı�Dirchilet�߽�������Ӧ�ĵ�
    bdidx = zeros(Ndof,1); 
    bdidx(fixedNode) = 1;
    Tbd = spdiags(bdidx,0,Ndof,Ndof);
    T = spdiags(1-bdidx,0,Ndof,Ndof);
    AD = T*A*T + Tbd;
    


    %2.��Neumann�߽�,�������Ҷ���b
    Neumann = bdEdge;
    el = sqrt(sum((node(bdEdge(:,1),:)-node(bdEdge(:,2),:)).^2,2));
    option.gNquadorder = 2;%%�������Ի���
    [lambdagN,weightgN] = quadpts1(option.gNquadorder);
    phigN = lambdagN;%�߽����Ի�����
    nQuadgN = size(lambdagN,1);
    ge = zeros(size(Neumann,1),2);
    for ppp = 1:nQuadgN
            ppxy = lambdagN(ppp,1)*node(Neumann(:,1),:) ...
                 + lambdagN(ppp,2)*node(Neumann(:,2),:);
            gNp = pde.g_N(ppxy);
            for igN = 1:2
                ge(:,igN) = ge(:,igN) + weightgN(ppp)*phigN(ppp,igN)*gNp;
            end
    end
    ge = ge.*repmat(el,1,2);
    b = b + accumarray(Neumann(:), ge(:),[Ndof,1]); 


    %%����Dirichlet�߽�������b;
    u = zeros(Ndof,1);
    u(fixedNode) = pde.g_D(node(fixedNode,:));
    b = b - A*u;
    b(fixedNode) = u(fixedNode);
    end
%%�����߽�����
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% HermiteP3
function [soln,eqn,info,area] = PoissonHermiteP3s(node,elem,bdFlag,pde,option)
%%2019������Ԫ�������ϻ���ҵ1,ΪP3HermiteԪ����[0,1]^2�ϵĻ�ϱ߽�������Possion����.
Nd = size(node,1); %�����;
NT = size(elem,1);%��Ԫ����;
[edge,elem2edge,elemdof] = dofP3Hermite(elem); %�ֲ���������Ӧ�������������λ��.
                                                          %��elemdof(i,j)����ʾ��i����Ԫ�ĵ�j����������Ӧ��������������; 
                                                          %edge(i,:)Ϊ��i���ߵ������˵����,elem2edge(i,j)Ϊ��i����Ԫ�ĵ�j����
                                                          %��ȫ�����е����.
                                                          
 Ndof = max(double(elemdof(:)));
 %%
%%��װ�նȾ���.
%%�ֲ�������������������Ķ�Ӧ��Ź�ϵ
time = cputime;
ii = zeros(55*NT,1);
jj = zeros(55*NT,1);
index = 0;
for i = 1:10
    for j = i:10
        ii(index+1:index+NT) = double(elemdof(:,i));
        jj(index+1:index+NT) = double(elemdof(:,j));
        index = index + NT;
    end
end
        
%Dlambda(i,:,s)��ʾ��i����Ԫ��lambda_s���ݶȵ���
[Dlambda,area] = gradbasis(node,elem);
option.fquadorder = 6;
[lambda,w] = quadpts(option.fquadorder);
nQuad = size(lambda,1);
sA = zeros(55*NT,nQuad);
%����ֲ��նȾ���
for p = 1:nQuad
    la = lambda(p,:);
    %�ֲ��������ݶ�
    Dphip = getDphis(la,Dlambda,node,elem);
    index = 0;
    for i = 1:10
        for j = i:10
            Aij = 0;    
            Aij =  Aij + w(p)*dot(Dphip(:,:,i),Dphip(:,:,j),2);
            Aij = Aij.*area;  
            sA(index+1:index+NT,p) =  Aij;
            index = index+NT;
        end
    end
end
sA = sum(sA,2);
diagIdx = (ii == jj);   
upperIdx = ~diagIdx;
A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),Ndof,Ndof);
AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),Ndof,Ndof);
A = A + AU + AU';
clear Aij ii jj sA




%%
%%�齨�Ҷ���b
[lambda,w] = quadpts(option.fquadorder); %4�׾��Ȼ���
nQuad = size(lambda,1);
bt = zeros(NT,10);
for p = 1:nQuad
    %����ֲ�������ֵ
     phi = getphis(lambda(p,:),node,elem);
     %���ֵ�
     pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        fp = pde.f(pxy);
        for j = 1:10
            bt(:,j) = bt(:,j) + w(p).*phi(:,j).*fp;
        end
end
bt = bt.*repmat(area,1,10);
b = accumarray(elemdof(:),bt(:),[Ndof,1]);

%% �߽���������
[AD,b,u,freeDof] = getbd(b);
%min(eig(full(AD)))
%%���u
%% Record assembeling time
assembleTime = cputime - time;
        tic;
        u(freeDof) = AD(freeDof,freeDof)\b(freeDof);
        residual = norm(b - AD*u);
        info = struct('solverTime',toc,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
%% Compute Du
Du = [];

%% Output
if nargout == 1
    soln = u;
else
    soln = struct('u',u,'Du',Du);
    eqn = struct('A',AD,'b',b,'edge',edge,'freeDof',freeDof,'Lap',A,'elemdof',elemdof);
    info.assembleTime = assembleTime;
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�Ӻ��� getbd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AD,b,u,freeDof] = getbd(b)
%%����Dirichle�߽�������Neumann�߽�����

%1. ��Dirichlet�߽�,����Ϊ:fixedNode,����������A
isFixedDof = false(Ndof,1);   
isFixedDof(elemdof(bdFlag(:,1)==1,[2,3,8,9])) = true;
isFixedDof(elemdof(bdFlag(:,2)==1,[1,3,7,9])) = true;
isFixedDof(elemdof(bdFlag(:,3)==1,[1,2,7,8])) = true;
fixedDof = find(isFixedDof);
freeDof = find(~isFixedDof);  






% Modify the matrix
% Build Dirichlet boundary condition into the matrix AD by enforcing
% AD(fixedDof,fixedDof)=I, AD(fixedDof,ufreeDof)=0, AD(ufreeDof,fixedDof)=0..
bdidx = zeros(Ndof,1); 
bdidx(fixedDof) = 1;
Tbd = spdiags(bdidx,0,Ndof,Ndof);
T = spdiags(1-bdidx,0,Ndof,Ndof);
AD = T*A*T + Tbd;


%%
%2.��Neumann�߽�
idxN = (bdFlag(:) == 2);      % all Neumann edges in bdFlag        
Neumannidx = elem2edge(idxN); % index of Neumann edges
% since boundary integral is also needed 
Neumann   = edge(Neumannidx,:);

%�Ľ׻���
 [lambdagN,weightgN] = quadpts1(option.fquadorder);
 nQuadgN = size(lambdagN,1);
 %�߽������on
 %(1--2),�ֱ�Ϊ������,����ƫ_{x}(.)(a_1),ƫ_{x}(.)(a_2),ƫ_{y}(.)(a_1),ƫ_{y}(.)(a_2);
el = sqrt(sum((node(Neumann(:,1),:) - node(Neumann(:,2),:)).^2,2));
ge = zeros(size(Neumann,1),6);
 
 for pp = 1:nQuadgN
         bdphi = getNeumannphi(lambdagN(pp,:),Neumann,node);
            ppxy = lambdagN(pp,1)*node(Neumann(:,1),:) ...
         + lambdagN(pp,2)*node(Neumann(:,2),:);
         gNp = pde.g_N(ppxy);
       
         for is = 1:6
            ge(:,is) = ge(:,is) + weightgN(pp)*gNp.*bdphi(:,is);
         end

 end

 
 
%%
%3 �����Ҷ��� 
%Neumann ����
ge = ge.*repmat(el,1,6);
b(1:3*Nd) = b(1:3*Nd) + accumarray([Neumann(:);Neumann(:)+Nd;Neumann(:)+2*Nd], ge(:),[3*Nd,1]);


% Dirichlet ����
u = zeros(Ndof,1);
idx = (fixedDof > Nd);                         % index of ƫ����
fixedNode = fixedDof(~idx);
u(fixedNode) = pde.g_D(node(fixedNode,:)); % bd value at vertex dofs

 Dg_d = pde.Dg_d(node(fixedNode,:));
 
 u(fixedDof(idx)) = Dg_d(:,2);

b = b - A*u;
b(fixedDof) = u(fixedDof);    




end


function [edge,elem2edge,elemdof] = dofP3Hermite(elem)
%%�ڸ�����Ԫ�ʷ�elem��������ֲ�������������������Ķ�Ӧ��ϵ
%%elemdof(i,k)��ʾ��i����Ԫ�ĵ�k���ֲ���������Ӧ��ȫ�ֻ����������.


NT = size(elem,1);%%��Ԫ�ʷָ���.
Nd = max(elem(:));%%�����


totalEdge = uint32([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])]);
totalEdges = sort(totalEdge,2);
[edge, ~, j] = myunique(totalEdges);
elem2edge = reshape(j,NT,3);


elemdof = [elem elem+Nd elem+2*Nd 3*Nd+(1:NT)'];

end

function Dphip = getDphis(la,Dla,node,elem)
a = zeros(size(elem,1),2,3);
a(:,:,1) = node(elem(:,1),:);
a(:,:,2) = node(elem(:,2),:);
a(:,:,3) = node(elem(:,3),:);

Dphip = zeros(size(elem,1),2,10);

Dla1la2la3 = la(1)*la(2).*Dla(:,:,3)+la(1).*Dla(:,:,2).*la(3)+Dla(:,:,1).*la(2).*la(3);


i = 1;
Dphip(:,:,1) = 3*la(i)^2.*Dla(:,:,i) + 6*la(i)*(la(i+1)+la(3)).*Dla(:,:,i)...
                            + 3*la(i)^2.*(Dla(:,:,i+1)+Dla(:,:,3))-7*Dla1la2la3;
                           
                           
                           
i = 2;
Dphip(:,:,2) = 3*la(i)^2.*Dla(:,:,i) + 6*la(i)*(la(i+1)+la(i-1)).*Dla(:,:,i)...
                            + 3*la(i)^2.*(Dla(:,:,i+1)+Dla(:,:,i-1))-7*Dla1la2la3;
                           
                          
                          
                          
i = 3;
Dphip(:,:,3) = 3*la(i)^2.*Dla(:,:,i) + 6*la(i)*(la(1)+la(i-1)).*Dla(:,:,i)...
                            + 3*la(i)^2.*(Dla(:,:,1)+Dla(:,:,i-1))-7*Dla1la2la3;
                          
                          
                          
i = 1;
idx = 4;


Dphip(:,:,idx) = 2*la(i).*(la(i+1).*(a(:,1,i+1)-a(:,1,i))+la(3).*(a(:,1,3)-a(:,1,i))) ...
                                  .*Dla(:,:,i)...
                                  +la(i)^2.*(Dla(:,:,i+1).*(a(:,1,i+1)-a(:,1,i))+Dla(:,:,3).*(a(:,1,3)-a(:,1,i)))...
                                  -(a(:,1,i+1)+a(:,1,3)-2*a(:,1,i)).*Dla1la2la3;
                              
Dphip(:,:,idx+3) = 2*la(i).*(la(i+1).*(a(:,2,i+1)-a(:,2,i))+la(3).*(a(:,2,3)-a(:,2,i))) ...
                                  .*Dla(:,:,i)...
                                  +la(i)^2.*(Dla(:,:,i+1).*(a(:,2,i+1)-a(:,2,i))+Dla(:,:,3).*(a(:,2,3)-a(:,2,i)))...
                                  -(a(:,2,i+1)+a(:,2,3)-2*a(:,2,i)).*Dla1la2la3; 
                              
                              
i = 2;
idx = 5;

  Dphip(:,:,idx) = 2*la(i).*(la(i+1).*(a(:,1,i+1)-a(:,1,i))+la(i-1).*(a(:,1,i-1)-a(:,1,i))) ...
                                  .*Dla(:,:,i)...
                                  +la(i)^2.*(Dla(:,:,i+1).*(a(:,1,i+1)-a(:,1,i))+Dla(:,:,i-1).*(a(:,1,i-1)-a(:,1,i)))...
                                  -(a(:,1,i+1)+a(:,1,i-1)-2*a(:,1,i)).*Dla1la2la3;
                              
Dphip(:,:,idx+3) = 2*la(i).*(la(i+1).*(a(:,2,i+1)-a(:,2,i))+la(i-1).*(a(:,2,i-1)-a(:,2,i))) ...
                                  .*Dla(:,:,i)...
                                  +la(i)^2.*(Dla(:,:,i+1).*(a(:,2,i+1)-a(:,2,i))+Dla(:,:,i-1).*(a(:,2,i-1)-a(:,2,i)))...
                                  -(a(:,2,i+1)+a(:,2,i-1)-2*a(:,2,i)).*Dla1la2la3;                               
                                    

                              

                              
                              
i = 3;
idx = 6;

  Dphip(:,:,idx) = 2*la(i).*(la(1).*(a(:,1,1)-a(:,1,i))+la(i-1).*(a(:,1,i-1)-a(:,1,i))) ...
                                  .*Dla(:,:,i)...
                                  +la(i)^2.*(Dla(:,:,1).*(a(:,1,1)-a(:,1,i))+Dla(:,:,i-1).*(a(:,1,i-1)-a(:,1,i)))...
                                  -(a(:,1,1)+a(:,1,i-1)-2*a(:,1,i)).*Dla1la2la3;
                              
Dphip(:,:,idx+3) = 2*la(i).*(la(1).*(a(:,2,1)-a(:,2,i))+la(i-1).*(a(:,2,i-1)-a(:,2,i))) ...
                                  .*Dla(:,:,i)...
                                  +la(i)^2.*(Dla(:,:,1).*(a(:,2,1)-a(:,2,i))+Dla(:,:,i-1).*(a(:,2,i-1)-a(:,2,i)))...
                                  -(a(:,2,1)+a(:,2,i-1)-2*a(:,2,i)).*Dla1la2la3;  

Dphip(:,:,10) = 27*Dla1la2la3;

end

function phi = getphis(la,node,elem)
a = zeros(size(elem,1),2,3);
a(:,:,1) = node(elem(:,1),:);
a(:,:,2) = node(elem(:,2),:);
a(:,:,3) = node(elem(:,3),:);



phi = zeros(size(elem,1),10);
la1la2la3 = la(1).*la(2).*la(3);


i = 1;
phi(:,1) = la(i).^3+3*la(i).^2.*(la(i+1)+la(3))-7*la1la2la3;
          
                     
i = 2;
phi(:,2) = la(i).^3+3*la(i).^2.*(la(i+1)+la(i-1))-7*la1la2la3;
                     
i = 3;
phi(:,3) = la(i).^3+3*la(i).^2.*(la(1)+la(i-1))-7*la1la2la3;


i = 1;
idx = 4;

phi(:,[idx,idx+3]) = la(i).^2.*(la(i+1).*(a(:,:,i+1)-a(:,:,i))+la(3).*(a(:,:,3)-a(:,:,i)))...
                                            -(a(:,:,i+1)+a(:,:,3)-2*a(:,:,i))*la1la2la3;
                                        
                                        
                                        
i = 2;
idx = 5;
phi(:,[idx,idx+3]) = la(i).^2.*(la(i+1).*(a(:,:,i+1)-a(:,:,i))+la(i-1).*(a(:,:,i-1)-a(:,:,i)))...
                                            -(a(:,:,i+1)+a(:,:,i-1)-2*a(:,:,i))*la1la2la3;
                                        
                                        
i = 3;
idx = 6;
phi(:,[idx,idx+3]) = la(i).^2.*(la(1).*(a(:,:,1)-a(:,:,i))+la(i-1).*(a(:,:,i-1)-a(:,:,i)))...
                                            -(a(:,:,1)+a(:,:,i-1)-2*a(:,:,i))*la1la2la3;
                                        
                                        
phi(:,10) = 27*la1la2la3;


% % % phi(:,idx) = 2.*area.*(la(i).^2.*(la(i+1).*e3(:,1)-la(3).*e2(:,1))...
% % %                                        +(e2(:,1)-e3(:,1)).*la1la2la3)./e_T23;
% % % phi(:,idx+3) = 2.*area.*(la(i).^2.*(la(i+1).*e3(:,2)-la(3).*e2(:,2))...
% % %                                        +(e2(:,2)-e3(:,2)).*la1la2la3)./e_T23;






end

function  bdphi = getNeumannphi(la,Neumann,node)
   a1 = node(Neumann(:,1),:);
   a2 = node(Neumann(:,2),:);
    bdphi = zeros(size(Neumann,1),6);
     i = 1;
     j = 2;
      bdphi(:,1) = la(i).^3+3*la(i).^2.*la(j);
      
      i = 2;
      j = 1;
       bdphi(:,2) =  la(i).^3+3*la(i).^2.*la(j);
       
       idx = 3;
       i = 1;
       j = 2;
       bdphi(:,[idx,idx+2]) = la(i)^2*la(j).*(a2-a1);
       
       idx = 4;
       i = 2;
       j = 1;
       bdphi(:,[idx,idx+2]) = la(i)^2*la(j).*(a1-a2);
       
       
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% HermiteP4
function [soln,eqn,info,area] = PoissonHermiteP4s(node,elem,bdFlag,pde,option)
%%2019������Ԫ�������ϻ���ҵ1,ΪP4HermiteԪ����[0,1]^2�ϵĻ�ϱ߽�������Possion����.
Nd = size(node,1); %�����;
NT = size(elem,1);%��Ԫ����;
[edge,elem2edge,elemdof] = dofP4Hermite(elem); %�ֲ���������Ӧ�������������λ��.
                                                          %��elemdof(i,j)����ʾ��i����Ԫ�ĵ�j����������Ӧ��������������; 
                                                          %edge(i,:)Ϊ��i���ߵ������˵����,elem2edge(i,j)Ϊ��i����Ԫ�ĵ�j����
                                                          %��ȫ�����е����.
 NE = size(edge,1);                                                         
 Ndof = max(double(elemdof(:)));
 %%
%%��װ�նȾ���.
%%�ֲ�������������������Ķ�Ӧ��Ź�ϵ
time = cputime;
ii = zeros(120*NT,1);
jj = zeros(120*NT,1);
index = 0;
for i = 1:15
    for j = i:15
        ii(index+1:index+NT) = double(elemdof(:,i));
        jj(index+1:index+NT) = double(elemdof(:,j));
        index = index + NT;
    end
end
        
%Dlambda(i,:,s)��ʾ��i����Ԫ��lambda_s���ݶȵ���
[Dlambda,area] = gradbasis(node,elem);
option.fquadorder = 8;
[lambda,w] = quadpts(option.fquadorder);
nQuad = size(lambda,1);
sA = zeros(120*NT,nQuad);
%����ֲ��նȾ���
for p = 1:nQuad
    la = lambda(p,:);
    %�ֲ��������ݶ�
    Dphip = getDphis(la,Dlambda,node,elem);
    index = 0;
    for i = 1:15
        for j = i:15
            Aij = 0;    
            Aij =  Aij + w(p)*dot(Dphip(:,:,i),Dphip(:,:,j),2);
            Aij = Aij.*area;  
            sA(index+1:index+NT,p) =  Aij;
            index = index+NT;
        end
    end
end
sA = sum(sA,2);
diagIdx = (ii == jj);   
upperIdx = ~diagIdx;
A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),Ndof,Ndof);
AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),Ndof,Ndof);
A = A + AU + AU';
clear Aij ii jj sA




%%
%%�齨�Ҷ���b
[lambda,w] = quadpts(option.fquadorder); 
nQuad = size(lambda,1);
bt = zeros(NT,15);
for p = 1:nQuad
    %����ֲ�������ֵ
     phi = getphis(lambda(p,:),node,elem);
     %���ֵ�
     pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:);
        fp = pde.f(pxy);
        for j = 1:15
            bt(:,j) = bt(:,j) + w(p).*phi(:,j).*fp;
        end
end
bt = bt.*repmat(area,1,15);
b = accumarray(elemdof(:),bt(:),[Ndof,1]);

%% �߽���������
[AD,b,u,freeDof] = getbd(b);
%min(eig(full(A)))
%%���u
%% Record assembeling time
assembleTime = cputime - time;
        tic;
        u(freeDof) = AD(freeDof,freeDof)\b(freeDof);
        residual = norm(b - AD*u);
        info = struct('solverTime',toc,'itStep',0,'err',residual,'flag',2,'stopErr',residual);
%% Compute Du
Du = [];

%% Output
if nargout == 1
    soln = u;
else
    soln = struct('u',u,'Du',Du);
    eqn = struct('A',AD,'b',b,'edge',edge,'freeDof',freeDof,'Lap',A,'elemdof',elemdof);
    info.assembleTime = assembleTime;
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�Ӻ��� getbd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [AD,b,u,freeDof] = getbd(b)
%%����Dirichle�߽�������Neumann�߽�����

%1. ��Dirichlet�߽�,����Ϊ:fixedNode,����������A
isFixedDof = false(Ndof,1);   
isFixedDof(elemdof(bdFlag(:,1)==1,[2,3,5,6,8,9,10])) = true;
isFixedDof(elemdof(bdFlag(:,2)==1,[1,3,4,6,7,9,11])) = true;
isFixedDof(elemdof(bdFlag(:,3)==1,[1,2,4,5,7,8,12])) = true;
fixedDof = find(isFixedDof);
freeDof = find(~isFixedDof);  






% Modify the matrix
% Build Dirichlet boundary condition into the matrix AD by enforcing
% AD(fixedDof,fixedDof)=I, AD(fixedDof,ufreeDof)=0, AD(ufreeDof,fixedDof)=0..
bdidx = zeros(Ndof,1); 
bdidx(fixedDof) = 1;
Tbd = spdiags(bdidx,0,Ndof,Ndof);
T = spdiags(1-bdidx,0,Ndof,Ndof);
AD = T*A*T + Tbd;


%%
%2.��Neumann�߽�
idxN = (bdFlag(:) == 2);      % all Neumann edges in bdFlag        
Neumannidx = elem2edge(idxN); % index of Neumann edges
% since boundary integral is also needed 
Neumann   = edge(Neumannidx,:);

%����
 [lambdagN,weightgN] = quadpts1(option.fquadorder);
 nQuadgN = size(lambdagN,1);
 %�߽������on
 %(1--2),�ֱ�Ϊ������,����ƫ_{x}(.)(a_1),ƫ_{x}(.)(a_2),ƫ_{y}(.)(a_1),ƫ_{y}(.)(a_2),�����е�ֵ
el = sqrt(sum((node(Neumann(:,1),:) - node(Neumann(:,2),:)).^2,2));
ge = zeros(size(Neumann,1),7);
 
 for pp = 1:nQuadgN
         bdphi = getNeumannphi(lambdagN(pp,:),Neumann,node);
            ppxy = lambdagN(pp,1)*node(Neumann(:,1),:) ...
         + lambdagN(pp,2)*node(Neumann(:,2),:);
         gNp = pde.g_N(ppxy);
       
         for is = 1:7
            ge(:,is) = ge(:,is) + weightgN(pp)*gNp.*bdphi(:,is);
         end

 end

 
 
%%
%3 �����Ҷ��� 
%Neumann ����
ge = ge.*repmat(el,1,7);
b(1:3*Nd+NE) = b(1:3*Nd+NE) + accumarray([Neumann(:);Neumann(:)+Nd;Neumann(:)+2*Nd;...
                               Neumannidx+3*Nd], ge(:),[3*Nd+NE,1]);


% Dirichlet ����
u = zeros(Ndof,1);
idx = (fixedDof <= Nd);                         % index of ��
fixedNode = fixedDof(idx);
u(fixedNode) = pde.g_D(node(fixedNode,:)); % bd value at vertex dofs

%%��ƫ����ֵ
Dg_d = pde.Dg_d(node(fixedNode,:));
idx = and(fixedDof<=3*Nd,fixedDof>Nd);

u(fixedDof(idx)) = Dg_d(:);




%%���е�ֵ
idx = (fixedDof>3*Nd);
bdEdgeIdx = fixedDof(idx) - 3*Nd;
bdEdgeMid = (node(edge(bdEdgeIdx,1),:) + node(edge(bdEdgeIdx,2),:))/2;
u(fixedDof(idx)) = pde.g_D(bdEdgeMid);



b = b - A*u;
b(fixedDof) = u(fixedDof);    




end





function [edge,elem2edge,elemdof] = dofP4Hermite(elem)
%%�ڸ�����Ԫ�ʷ�elem��������ֲ�������������������Ķ�Ӧ��ϵ
%%elemdof(i,k)��ʾ��i����Ԫ�ĵ�k���ֲ���������Ӧ��ȫ�ֻ����������.


NT = size(elem,1);%%��Ԫ�ʷָ���.
Nd = max(elem(:));%%�����


totalEdge = uint32([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])]);
totalEdges = sort(totalEdge,2);
[edge, ~, j] = myunique(totalEdges);
elem2edge = reshape(j,NT,3);

NE = size(edge,1);
elemdof = [elem elem+Nd elem+2*Nd 3*Nd+elem2edge ...
                          3*Nd+NE+(1:NT)' 3*Nd+NE+NT+(1:NT)' 3*Nd+NE+2*NT+(1:NT)'];

end

function Dphip = getDphis(la,Dla,node,elem)
a = zeros(size(elem,1),2,3);
a(:,:,1) = node(elem(:,1),:);
a(:,:,2) = node(elem(:,2),:);
a(:,:,3) = node(elem(:,3),:);

Dphip = zeros(size(elem,1),2,15);
la1la2la3 = la(1)*la(2)*la(3);
Dla1la2la3 = la(1)*la(2)*Dla(:,:,3)+la(1)*Dla(:,:,2)*la(3)+Dla(:,:,1)*la(2)*la(3);
suma = la(1).*a(:,:,1)+la(2).*a(:,:,2)+la(3).*a(:,:,3);
Dsuma1 =a(:,1,1).*Dla(:,:,1)+a(:,1,2).*Dla(:,:,2)+a(:,1,3).*Dla(:,:,3);
Dsuma2 =a(:,2,1).*Dla(:,:,1)+a(:,2,2).*Dla(:,:,2)+a(:,2,3).*Dla(:,:,3);

i = 1;
Dphip(:,:,1) = 12*la(i)^2*(1-la(i)).*Dla(:,:,i)...
                                 -10*la(i)*(la(3)^2+la(i+1)^2).*Dla(:,:,i)...
                                 -10*la(i)^2.*(la(3).*Dla(:,:,3)+la(i+1).*Dla(:,:,i+1))...
                                 +(8-26*la(i)).*Dla1la2la3 -la1la2la3.*26*Dla(:,:,i);
                            
                           
i = 2;
Dphip(:,:,2) = 12*la(i)^2*(1-la(i)).*Dla(:,:,i)...
                                 -10*la(i)*(la(i-1)^2+la(i+1)^2).*Dla(:,:,i)...
                                 -10*la(i)^2.*(la(i-1).*Dla(:,:,i-1)+la(i+1).*Dla(:,:,i+1))...
                                 +(8-26*la(i)).*Dla1la2la3+la1la2la3*(-26*Dla(:,:,i));
                            
                           
                          
                          
                          
i = 3;
Dphip(:,:,3) = 12*la(i)^2*(1-la(i)).*Dla(:,:,i)...
                                 -10*la(i)*(la(i-1)^2+la(1)^2).*Dla(:,:,i)...
                                 -10*la(i)^2.*(la(i-1).*Dla(:,:,i-1)+la(1).*Dla(:,:,1))...
                                 +(8-26*la(i)).*Dla1la2la3+la1la2la3*(-26*Dla(:,:,i));
                          

                             
i = 1;
Dphip(:,:,10) = Dla1la2la3.*(24*la(i)-14)+24*la1la2la3.*Dla(:,:,i)...
                                +32*(la(3)*la(i+1)^2.*Dla(:,:,3)+la(i+1)*la(3)^2.*Dla(:,:,i+1));
                            
                            
                            
i = 2;
Dphip(:,:,11) = Dla1la2la3.*(24*la(i)-14)+24*la1la2la3.*Dla(:,:,i)...
                                +32*(la(i-1)*la(i+1)^2.*Dla(:,:,i-1)+la(i+1)*la(i-1)^2.*Dla(:,:,i+1));


i = 3;                           
Dphip(:,:,12) = Dla1la2la3.*(24*la(i)-14)+24*la1la2la3.*Dla(:,:,i)...
                                +32*(la(i-1)*la(1)^2.*Dla(:,:,i-1)+la(1)*la(i-1)^2.*Dla(:,:,1));
                            
                            
                            
                            
for i = 1:3
    Dphip(:,:,12+i)  = 32*Dla1la2la3.*(4*la(i)-1)+128*la1la2la3.*Dla(:,:,i);
end
                            
                            
                            
                            
                            
                            
i = 1;
idx = 4;


Dphip(:,:,idx) = 3*la(i)^2*(suma(:,1)-a(:,1,i)).*Dla(:,:,i)+la(i)^3.*Dsuma1 ...
                                 +2*la(i)*(a(:,1,i).*(la(3)^2+la(i+1)^2)-(a(:,1,3)*la(3)^2 ...
                                  + a(:,1,i+1)*la(i+1)^2)).*Dla(:,:,i)...
                                  +2*la(i)^2.*(a(:,1,i).*(la(3).*Dla(:,:,3)+la(i+1).*Dla(:,:,i+1))...
                                  - (a(:,1,3).*la(3).*Dla(:,:,3)+a(:,1,i+1).*la(i+1).*Dla(:,:,i+1)))...
                                  +0.25*Dla1la2la3.*(a(:,1,i).*(20*la(i)-6)+a(:,1,3).*(4*la(3)...
                                  -8*la(i)+1)+a(:,1,i+1).*(4*la(i+1)-8*la(i)+1))...
                                  +0.25*la1la2la3.*(20*a(:,1,i).*Dla(:,:,i)...
                                  +a(:,1,3).*(4*Dla(:,:,3)-8*Dla(:,:,i))...
                                  +a(:,1,i+1).*(4*Dla(:,:,i+1)-8*Dla(:,:,i)));
                              
Dphip(:,:,idx+3) = 3*la(i)^2*(suma(:,2)-a(:,2,i)).*Dla(:,:,i)+la(i)^3.*Dsuma2 ...
                                 +2*la(i)*(a(:,2,i).*(la(3)^2+la(i+1)^2)-(a(:,2,3)*la(3)^2 ...
                                  + a(:,2,i+1)*la(i+1)^2)).*Dla(:,:,i)...
                                  +2*la(i)^2.*(a(:,2,i).*(la(3).*Dla(:,:,3)+la(i+1).*Dla(:,:,i+1))...
                                  - (a(:,2,3).*la(3).*Dla(:,:,3)+a(:,2,i+1).*la(i+1).*Dla(:,:,i+1)))...
                                  +0.25*Dla1la2la3.*(a(:,2,i).*(20*la(i)-6)+a(:,2,3).*(4*la(3)...
                                  -8*la(i)+1)+a(:,2,i+1).*(4*la(i+1)-8*la(i)+1))...
                                  +0.25*la1la2la3.*(20*a(:,2,i).*Dla(:,:,i)...
                                  +a(:,2,3).*(4*Dla(:,:,3)-8*Dla(:,:,i))...
                                  +a(:,2,i+1).*(4*Dla(:,:,i+1)-8*Dla(:,:,i)));
                              
i = 2;
idx = 5;


Dphip(:,:,idx) = 3*la(i)^2*(suma(:,1)-a(:,1,i)).*Dla(:,:,i)+la(i)^3.*Dsuma1 ...
                                 +2*la(i)*(a(:,1,i).*(la(i-1)^2+la(i+1)^2)-(a(:,1,i-1)*la(i-1)^2 ...
                                  + a(:,1,i+1)*la(i+1)^2)).*Dla(:,:,i)...
                                  +2*la(i)^2.*(a(:,1,i).*(la(i-1).*Dla(:,:,i-1)+la(i+1).*Dla(:,:,i+1))...
                                  - (a(:,1,i-1).*la(i-1).*Dla(:,:,i-1)+a(:,1,i+1).*la(i+1).*Dla(:,:,i+1)))...
                                  +0.25*Dla1la2la3.*(a(:,1,i).*(20*la(i)-6)+a(:,1,i-1).*(4*la(i-1)...
                                  -8*la(i)+1)+a(:,1,i+1).*(4*la(i+1)-8*la(i)+1))...
                                  +0.25*la1la2la3.*(20*a(:,1,i).*Dla(:,:,i)...
                                  +a(:,1,i-1).*(4*Dla(:,:,i-1)-8*Dla(:,:,i))...
                                  +a(:,1,i+1).*(4*Dla(:,:,i+1)-8*Dla(:,:,i)));
                              
Dphip(:,:,idx+3) = 3*la(i)^2*(suma(:,2)-a(:,2,i)).*Dla(:,:,i)+la(i)^3.*Dsuma2 ...
                                 +2*la(i)*(a(:,2,i).*(la(i-1)^2+la(i+1)^2)-(a(:,2,i-1)*la(i-1)^2 ...
                                  + a(:,2,i+1)*la(i+1)^2)).*Dla(:,:,i)...
                                  +2*la(i)^2.*(a(:,2,i).*(la(i-1).*Dla(:,:,i-1)+la(i+1).*Dla(:,:,i+1))...
                                  - (a(:,2,i-1).*la(i-1).*Dla(:,:,i-1)+a(:,2,i+1).*la(i+1).*Dla(:,:,i+1)))...
                                  +0.25*Dla1la2la3.*(a(:,2,i).*(20*la(i)-6)+a(:,2,i-1).*(4*la(i-1)...
                                  -8*la(i)+1)+a(:,2,i+1).*(4*la(i+1)-8*la(i)+1))...
                                  +0.25*la1la2la3.*(20*a(:,2,i).*Dla(:,:,i)...
                                  +a(:,2,i-1).*(4*Dla(:,:,i-1)-8*Dla(:,:,i))...
                                  +a(:,2,i+1).*(4*Dla(:,:,i+1)-8*Dla(:,:,i)));
                              
                              
                              
                              
i = 3;
idx = 6;


Dphip(:,:,idx) = 3*la(i)^2*(suma(:,1)-a(:,1,i)).*Dla(:,:,i)+la(i)^3.*Dsuma1 ...
                                 +2*la(i)*(a(:,1,i).*(la(i-1)^2+la(1)^2)-(a(:,1,i-1)*la(i-1)^2 ...
                                  + a(:,1,1)*la(1)^2)).*Dla(:,:,i)...
                                  +2*la(i)^2.*(a(:,1,i).*(la(i-1).*Dla(:,:,i-1)+la(1).*Dla(:,:,1))...
                                  - (a(:,1,i-1).*la(i-1).*Dla(:,:,i-1)+a(:,1,1).*la(1).*Dla(:,:,1)))...
                                  +0.25*Dla1la2la3.*(a(:,1,i).*(20*la(i)-6)+a(:,1,i-1).*(4*la(i-1)...
                                  -8*la(i)+1)+a(:,1,1).*(4*la(1)-8*la(i)+1))...
                                  +0.25*la1la2la3.*(20*a(:,1,i).*Dla(:,:,i)...
                                  +a(:,1,i-1).*(4*Dla(:,:,i-1)-8*Dla(:,:,i))...
                                  +a(:,1,1).*(4*Dla(:,:,1)-8*Dla(:,:,i)));
                              
Dphip(:,:,idx+3) = 3*la(i)^2*(suma(:,2)-a(:,2,i)).*Dla(:,:,i)+la(i)^3.*Dsuma2 ...
                                 +2*la(i)*(a(:,2,i).*(la(i-1)^2+la(1)^2)-(a(:,2,i-1)*la(i-1)^2 ...
                                  + a(:,2,1)*la(1)^2)).*Dla(:,:,i)...
                                  +2*la(i)^2.*(a(:,2,i).*(la(i-1).*Dla(:,:,i-1)+la(1).*Dla(:,:,1))...
                                  - (a(:,2,i-1).*la(i-1).*Dla(:,:,i-1)+a(:,2,1).*la(1).*Dla(:,:,1)))...
                                  +0.25*Dla1la2la3.*(a(:,2,i).*(20*la(i)-6)+a(:,2,i-1).*(4*la(i-1)...
                                  -8*la(i)+1)+a(:,2,1).*(4*la(1)-8*la(i)+1))...
                                  +0.25*la1la2la3.*(20*a(:,2,i).*Dla(:,:,i)...
                                  +a(:,2,i-1).*(4*Dla(:,:,i-1)-8*Dla(:,:,i))...
                                  +a(:,2,1).*(4*Dla(:,:,1)-8*Dla(:,:,i)));
                              
                              
                              
                              
                              
                              
                              
end

function phi = getphis(la,node,elem)
a = zeros(size(elem,1),2,3);
a(:,:,1) = node(elem(:,1),:);
a(:,:,2) = node(elem(:,2),:);
a(:,:,3) = node(elem(:,3),:);



phi = zeros(size(elem,1),15);
la1la2la3 = la(1).*la(2).*la(3);
suma = la(1).*a(:,:,1)+la(2).*a(:,:,2)+la(3).*a(:,:,3);

i = 1;
phi(:,1) = 4*la(i).^3-3*la(i).^4-5*la(i).^2.*(la(i+1).^2+la(3).^2)+...
                     la1la2la3.*(8-26*la(i));
                 
          
                     
i = 2;
phi(:,2) = 4*la(i).^3-3*la(i).^4-5*la(i).^2.*(la(i+1).^2+la(i-1).^2)+...
                     la1la2la3.*(8-26*la(i));
                     
i = 3;
phi(:,3) = 4*la(i).^3-3*la(i).^4-5*la(i).^2.*(la(1).^2+la(i-1).^2)+...
                     la1la2la3.*(8-26*la(i));


                 
                 
                 
i = 1;
phi(:,10) = la1la2la3.*(24*la(i)-14)+16*la(i+1).^2.*la(3).^2;


i = 2;
phi(:,11) = la1la2la3.*(24*la(i)-14)+16*la(i+1).^2.*la(i-1).^2;


i = 3;
phi(:,12) = la1la2la3.*(24*la(i)-14)+16*la(1).^2.*la(i-1).^2;

                 
                 
for i = 1:3
  phi(:,12+i) =32* la1la2la3.*(4*la(i)-1);
end               
  


i = 1;
idx = 4;

phi(:,[idx,idx+3]) = la(i).^3.*(suma-a(:,:,i))+la(i).^2.*((la(i+1).^2 ...
                                          +la(3).^2).*a(:,:,i)-(la(i+1).^2.*a(:,:,i+1)+...
                                          la(3).^2.*a(:,:,3)))...
                                          +0.25*la1la2la3.*(a(:,:,i).*(20*la(i)-6)+a(:,:,3).*...
                                          (4*la(3)-8*la(i)+1)+a(:,:,i+1).*(4*la(i+1)-8*la(i)+1));
                                      
i = 2;
idx = 5;                                      
 
phi(:,[idx,idx+3]) = la(i).^3.*(suma-a(:,:,i))+la(i).^2.*((la(i+1).^2 ...
                                          +la(i-1).^2).*a(:,:,i)-(la(i+1).^2.*a(:,:,i+1)+...
                                          la(i-1).^2.*a(:,:,i-1)))...
                                          +la1la2la3.*(a(:,:,i).*(20*la(i)-6)+a(:,:,i-1).*...
                                          (4*la(i-1)-8*la(i)+1)+a(:,:,i+1).*(4*la(i+1)-8*la(i)+1))/4;
                                      
i = 3;
idx = 6;                                      
 
phi(:,[idx,idx+3]) = la(i).^3.*(suma-a(:,:,i))+la(i).^2.*((la(1).^2 ...
                                          +la(i-1).^2).*a(:,:,i)-(la(1).^2.*a(:,:,1)+...
                                          la(i-1).^2.*a(:,:,i-1)))...
                                          +la1la2la3.*(a(:,:,i).*(20*la(i)-6)+a(:,:,i-1).*...
                                          (4*la(i-1)-8*la(i)+1)+a(:,:,1).*(4*la(1)-8*la(i)+1))/4;
                                                                              


end

function  bdphi = getNeumannphi(la,Neumann,node)
   a1 = node(Neumann(:,1),:);
   a2 = node(Neumann(:,2),:);
    bdphi = zeros(size(Neumann,1),7);
     i = 1;
     j = 2;
      bdphi(:,1) = la(i).^3.*(4-3*la(i))-5*la(i).^2.*la(j).^2;
      
      i = 2;
      j = 1;
       bdphi(:,2) =  la(i).^3.*(4-3*la(i))-5*la(i).^2.*la(j).^2;
       
       idx = 3;
       i = 1;
       j = 2;
       bdphi(:,[idx,idx+2]) = la(i).^3.*(la(1).*a1+la(2).*a2-a1)+la(i).^2.*la(j).^2.*(a1-a2);
       
       idx = 4;
       i = 2;
       j = 1;
       bdphi(:,[idx,idx+2]) = la(i).^3.*(la(1).*a1+la(2).*a2-a2)+la(i).^2.*la(j).^2.*(a2-a1);
       
       
       bdphi(:,7) = 16.*la(1).^2.*la(2).^2;
       
end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% LagrangeP1 Linf���
function err = getLinferror(node,elem,exactu,uh)
%%�õ�P1��Ƭ����Ԫ��Linf���
lambda = [1/3,1/3,1/3];
phi = lambda;
uhp = uh(elem(:,1))*phi(:,1)+uh(elem(:,2))*phi(:,2)+...
                     uh(elem(:,3))*phi(:,3);

pxy = node(elem(:,1),:)*lambda(1)+node(elem(:,2),:)*lambda(2)+...
                       node(elem(:,3),:)*lambda(3);

up = exactu(pxy);

err = max(abs(up-uhp));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%HermiteP3 Linf���
function err = getLinfH3errors(node,elemdof,elem,exactu,uh)
%%�õ�Hermite P3Ԫ��Linf���
option.fquadorder = 6;
err = 0;
[lambda,w] = quadpts(option.fquadorder);
nQuad = size(lambda,1);
for p = 1:nQuad
        phi = getphis(lambda(p,:),node,elem);
    uhp = 0;

    for i = 1:10
        uhp = uhp +uh(elemdof(:,i)).*phi(:,i);
    end




    pxy = node(elem(:,1),:)*lambda(p,1)+node(elem(:,2),:)*lambda(p,2)+...
                       node(elem(:,3),:)*lambda(p,3);

    up = exactu(pxy);

    err = max([(abs(up-uhp));err]);
end





function phi = getphis(la,node,elem)
a = zeros(size(elem,1),2,3);
a(:,:,1) = node(elem(:,1),:);
a(:,:,2) = node(elem(:,2),:);
a(:,:,3) = node(elem(:,3),:);



phi = zeros(size(elem,1),10);
la1la2la3 = la(1).*la(2).*la(3);


i = 1;
phi(:,1) = la(i).^3+3*la(i).^2.*(la(i+1)+la(3))-7*la1la2la3;
          
                     
i = 2;
phi(:,2) = la(i).^3+3*la(i).^2.*(la(i+1)+la(i-1))-7*la1la2la3;
                     
i = 3;
phi(:,3) = la(i).^3+3*la(i).^2.*(la(1)+la(i-1))-7*la1la2la3;


i = 1;
idx = 4;

phi(:,[idx,idx+3]) = la(i).^2.*(la(i+1).*(a(:,:,i+1)-a(:,:,i))+la(3).*(a(:,:,3)-a(:,:,i)))...
                                            -(a(:,:,i+1)+a(:,:,3)-2*a(:,:,i))*la1la2la3;
                                        
                                        
                                        
i = 2;
idx = 5;
phi(:,[idx,idx+3]) = la(i).^2.*(la(i+1).*(a(:,:,i+1)-a(:,:,i))+la(i-1).*(a(:,:,i-1)-a(:,:,i)))...
                                            -(a(:,:,i+1)+a(:,:,i-1)-2*a(:,:,i))*la1la2la3;
                                        
                                        
i = 3;
idx = 6;
phi(:,[idx,idx+3]) = la(i).^2.*(la(1).*(a(:,:,1)-a(:,:,i))+la(i-1).*(a(:,:,i-1)-a(:,:,i)))...
                                            -(a(:,:,1)+a(:,:,i-1)-2*a(:,:,i))*la1la2la3;
                                        
                                        
phi(:,10) = 27*la1la2la3;


% % % phi(:,idx) = 2.*area.*(la(i).^2.*(la(i+1).*e3(:,1)-la(3).*e2(:,1))...
% % %                                        +(e2(:,1)-e3(:,1)).*la1la2la3)./e_T23;
% % % phi(:,idx+3) = 2.*area.*(la(i).^2.*(la(i+1).*e3(:,2)-la(3).*e2(:,2))...
% % %                                        +(e2(:,2)-e3(:,2)).*la1la2la3)./e_T23;






end




end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%HermiteP4 Linf���
function err = getLinfH4errors(node,elemdof,elem,exactu,uh)
%%�õ�Hermite P4Ԫ��Linf���
option.fquadorder = 6;
err = 0;
[lambda,w] = quadpts(option.fquadorder);
nQuad = size(lambda,1);
for p = 1:nQuad
   
        phi = getphis(lambda(p,:),node,elem);


    uhp = 0;

    for i = 1:15
        uhp = uhp +uh(elemdof(:,i)).*phi(:,i);
    end




    pxy = node(elem(:,1),:)*lambda(p,1)+node(elem(:,2),:)*lambda(p,2)+...
                       node(elem(:,3),:)*lambda(p,3);

    up = exactu(pxy);

    err = max([(abs(up-uhp));err]);
end







function phi = getphis(la,node,elem)
a = zeros(size(elem,1),2,3);
a(:,:,1) = node(elem(:,1),:);
a(:,:,2) = node(elem(:,2),:);
a(:,:,3) = node(elem(:,3),:);



phi = zeros(size(elem,1),15);
la1la2la3 = la(1).*la(2).*la(3);
suma = la(1).*a(:,:,1)+la(2).*a(:,:,2)+la(3).*a(:,:,3);

i = 1;
phi(:,1) = 4*la(i).^3-3*la(i).^4-5*la(i).^2.*(la(i+1).^2+la(3).^2)+...
                     la1la2la3.*(8-26*la(i));
                 
          
                     
i = 2;
phi(:,2) = 4*la(i).^3-3*la(i).^4-5*la(i).^2.*(la(i+1).^2+la(i-1).^2)+...
                     la1la2la3.*(8-26*la(i));
                     
i = 3;
phi(:,3) = 4*la(i).^3-3*la(i).^4-5*la(i).^2.*(la(1).^2+la(i-1).^2)+...
                     la1la2la3.*(8-26*la(i));


                 
                 
                 
i = 1;
phi(:,10) = la1la2la3.*(24*la(i)-14)+16*la(i+1).^2.*la(3).^2;


i = 2;
phi(:,11) = la1la2la3.*(24*la(i)-14)+16*la(i+1).^2.*la(i-1).^2;


i = 3;
phi(:,12) = la1la2la3.*(24*la(i)-14)+16*la(1).^2.*la(i-1).^2;

                 
                 
for i = 1:3
  phi(:,12+i) =32* la1la2la3.*(4*la(i)-1);
end               
  


i = 1;
idx = 4;

phi(:,[idx,idx+3]) = la(i).^3.*(suma-a(:,:,i))+la(i).^2.*((la(i+1).^2 ...
                                          +la(3).^2).*a(:,:,i)-(la(i+1).^2.*a(:,:,i+1)+...
                                          la(3).^2.*a(:,:,3)))...
                                          +0.25*la1la2la3.*(a(:,:,i).*(20*la(i)-6)+a(:,:,3).*...
                                          (4*la(3)-8*la(i)+1)+a(:,:,i+1).*(4*la(i+1)-8*la(i)+1));
                                      
i = 2;
idx = 5;                                      
 
phi(:,[idx,idx+3]) = la(i).^3.*(suma-a(:,:,i))+la(i).^2.*((la(i+1).^2 ...
                                          +la(i-1).^2).*a(:,:,i)-(la(i+1).^2.*a(:,:,i+1)+...
                                          la(i-1).^2.*a(:,:,i-1)))...
                                          +la1la2la3.*(a(:,:,i).*(20*la(i)-6)+a(:,:,i-1).*...
                                          (4*la(i-1)-8*la(i)+1)+a(:,:,i+1).*(4*la(i+1)-8*la(i)+1))/4;
                                      
i = 3;
idx = 6;                                      
 
phi(:,[idx,idx+3]) = la(i).^3.*(suma-a(:,:,i))+la(i).^2.*((la(1).^2 ...
                                          +la(i-1).^2).*a(:,:,i)-(la(1).^2.*a(:,:,1)+...
                                          la(i-1).^2.*a(:,:,i-1)))...
                                          +la1la2la3.*(a(:,:,i).*(20*la(i)-6)+a(:,:,i-1).*...
                                          (4*la(i-1)-8*la(i)+1)+a(:,:,1).*(4*la(1)-8*la(i)+1))/4;
                                                                              


end






end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%HermiteP3 H1���
function err = getH1H3errors(node,elemdof,elem,exactDu,uh)
%P3HermiteԪ��H1�뷶���
[Dla,area] = gradbasis(node,elem);
option.fquadorder = 6;
err = 0;
[lambda,w] = quadpts(option.fquadorder);
nQuad = size(lambda,1);
for p = 1:nQuad
        Dphip = getDphis(lambda(p,:),Dla,node,elem);
       Duhp = 0;

    for i = 1:10
        Duhp = Duhp +uh(elemdof(:,i)).*Dphip(:,:,i);
    end

    pxy = node(elem(:,1),:)*lambda(p,1)+node(elem(:,2),:)*lambda(p,2)+...
                       node(elem(:,3),:)*lambda(p,3);

    Dup = exactDu(pxy);

    err = err + w(p)*sum((Dup-Duhp).^2,2);
end

err = err.*area;
err = sqrt(sum(err));








function Dphip = getDphis(la,Dla,node,elem)
a = zeros(size(elem,1),2,3);
a(:,:,1) = node(elem(:,1),:);
a(:,:,2) = node(elem(:,2),:);
a(:,:,3) = node(elem(:,3),:);

Dphip = zeros(size(elem,1),2,10);

Dla1la2la3 = la(1)*la(2).*Dla(:,:,3)+la(1).*Dla(:,:,2).*la(3)+Dla(:,:,1).*la(2).*la(3);


i = 1;
Dphip(:,:,1) = 3*la(i)^2.*Dla(:,:,i) + 6*la(i)*(la(i+1)+la(3)).*Dla(:,:,i)...
                            + 3*la(i)^2.*(Dla(:,:,i+1)+Dla(:,:,3))-7*Dla1la2la3;
                           
                           
                           
i = 2;
Dphip(:,:,2) = 3*la(i)^2.*Dla(:,:,i) + 6*la(i)*(la(i+1)+la(i-1)).*Dla(:,:,i)...
                            + 3*la(i)^2.*(Dla(:,:,i+1)+Dla(:,:,i-1))-7*Dla1la2la3;
                           
                          
                          
                          
i = 3;
Dphip(:,:,3) = 3*la(i)^2.*Dla(:,:,i) + 6*la(i)*(la(1)+la(i-1)).*Dla(:,:,i)...
                            + 3*la(i)^2.*(Dla(:,:,1)+Dla(:,:,i-1))-7*Dla1la2la3;
                          
                          
                          
i = 1;
idx = 4;


Dphip(:,:,idx) = 2*la(i).*(la(i+1).*(a(:,1,i+1)-a(:,1,i))+la(3).*(a(:,1,3)-a(:,1,i))) ...
                                  .*Dla(:,:,i)...
                                  +la(i)^2.*(Dla(:,:,i+1).*(a(:,1,i+1)-a(:,1,i))+Dla(:,:,3).*(a(:,1,3)-a(:,1,i)))...
                                  -(a(:,1,i+1)+a(:,1,3)-2*a(:,1,i)).*Dla1la2la3;
                              
Dphip(:,:,idx+3) = 2*la(i).*(la(i+1).*(a(:,2,i+1)-a(:,2,i))+la(3).*(a(:,2,3)-a(:,2,i))) ...
                                  .*Dla(:,:,i)...
                                  +la(i)^2.*(Dla(:,:,i+1).*(a(:,2,i+1)-a(:,2,i))+Dla(:,:,3).*(a(:,2,3)-a(:,2,i)))...
                                  -(a(:,2,i+1)+a(:,2,3)-2*a(:,2,i)).*Dla1la2la3; 
                              
                              
i = 2;
idx = 5;

  Dphip(:,:,idx) = 2*la(i).*(la(i+1).*(a(:,1,i+1)-a(:,1,i))+la(i-1).*(a(:,1,i-1)-a(:,1,i))) ...
                                  .*Dla(:,:,i)...
                                  +la(i)^2.*(Dla(:,:,i+1).*(a(:,1,i+1)-a(:,1,i))+Dla(:,:,i-1).*(a(:,1,i-1)-a(:,1,i)))...
                                  -(a(:,1,i+1)+a(:,1,i-1)-2*a(:,1,i)).*Dla1la2la3;
                              
Dphip(:,:,idx+3) = 2*la(i).*(la(i+1).*(a(:,2,i+1)-a(:,2,i))+la(i-1).*(a(:,2,i-1)-a(:,2,i))) ...
                                  .*Dla(:,:,i)...
                                  +la(i)^2.*(Dla(:,:,i+1).*(a(:,2,i+1)-a(:,2,i))+Dla(:,:,i-1).*(a(:,2,i-1)-a(:,2,i)))...
                                  -(a(:,2,i+1)+a(:,2,i-1)-2*a(:,2,i)).*Dla1la2la3;                               
                                    

                              

                              
                              
i = 3;
idx = 6;

  Dphip(:,:,idx) = 2*la(i).*(la(1).*(a(:,1,1)-a(:,1,i))+la(i-1).*(a(:,1,i-1)-a(:,1,i))) ...
                                  .*Dla(:,:,i)...
                                  +la(i)^2.*(Dla(:,:,1).*(a(:,1,1)-a(:,1,i))+Dla(:,:,i-1).*(a(:,1,i-1)-a(:,1,i)))...
                                  -(a(:,1,1)+a(:,1,i-1)-2*a(:,1,i)).*Dla1la2la3;
                              
Dphip(:,:,idx+3) = 2*la(i).*(la(1).*(a(:,2,1)-a(:,2,i))+la(i-1).*(a(:,2,i-1)-a(:,2,i))) ...
                                  .*Dla(:,:,i)...
                                  +la(i)^2.*(Dla(:,:,1).*(a(:,2,1)-a(:,2,i))+Dla(:,:,i-1).*(a(:,2,i-1)-a(:,2,i)))...
                                  -(a(:,2,1)+a(:,2,i-1)-2*a(:,2,i)).*Dla1la2la3;  

Dphip(:,:,10) = 27*Dla1la2la3;

end




end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%HermiteP4 H1���
function err = getH1H4errors(node,elemdof,elem,exactDu,uh)
%P4HermiteԪ��H1�뷶���
[Dla,area] = gradbasis(node,elem);
option.fquadorder = 6;
err = 0;
[lambda,w] = quadpts(option.fquadorder);
nQuad = size(lambda,1);
for p = 1:nQuad
        Dphip = getDphis(lambda(p,:),Dla,node,elem);
       Duhp = 0;

    for i = 1:15
        Duhp = Duhp +uh(elemdof(:,i)).*Dphip(:,:,i);
    end

    pxy = node(elem(:,1),:)*lambda(p,1)+node(elem(:,2),:)*lambda(p,2)+...
                       node(elem(:,3),:)*lambda(p,3);

    Dup = exactDu(pxy);

    err = err + w(p)*sum((Dup-Duhp).^2,2);
end

err = err.*area;
err = sqrt(sum(err));







function Dphip = getDphis(la,Dla,node,elem)
a = zeros(size(elem,1),2,3);
a(:,:,1) = node(elem(:,1),:);
a(:,:,2) = node(elem(:,2),:);
a(:,:,3) = node(elem(:,3),:);

Dphip = zeros(size(elem,1),2,15);
la1la2la3 = la(1)*la(2)*la(3);
Dla1la2la3 = la(1)*la(2)*Dla(:,:,3)+la(1)*Dla(:,:,2)*la(3)+Dla(:,:,1)*la(2)*la(3);
suma = la(1).*a(:,:,1)+la(2).*a(:,:,2)+la(3).*a(:,:,3);
Dsuma1 =a(:,1,1).*Dla(:,:,1)+a(:,1,2).*Dla(:,:,2)+a(:,1,3).*Dla(:,:,3);
Dsuma2 =a(:,2,1).*Dla(:,:,1)+a(:,2,2).*Dla(:,:,2)+a(:,2,3).*Dla(:,:,3);

i = 1;
Dphip(:,:,1) = 12*la(i)^2*(1-la(i)).*Dla(:,:,i)...
                                 -10*la(i)*(la(3)^2+la(i+1)^2).*Dla(:,:,i)...
                                 -10*la(i)^2.*(la(3).*Dla(:,:,3)+la(i+1).*Dla(:,:,i+1))...
                                 +(8-26*la(i)).*Dla1la2la3 -la1la2la3.*26*Dla(:,:,i);
                            
                           
i = 2;
Dphip(:,:,2) = 12*la(i)^2*(1-la(i)).*Dla(:,:,i)...
                                 -10*la(i)*(la(i-1)^2+la(i+1)^2).*Dla(:,:,i)...
                                 -10*la(i)^2.*(la(i-1).*Dla(:,:,i-1)+la(i+1).*Dla(:,:,i+1))...
                                 +(8-26*la(i)).*Dla1la2la3+la1la2la3*(-26*Dla(:,:,i));
                            
                           
                          
                          
                          
i = 3;
Dphip(:,:,3) = 12*la(i)^2*(1-la(i)).*Dla(:,:,i)...
                                 -10*la(i)*(la(i-1)^2+la(1)^2).*Dla(:,:,i)...
                                 -10*la(i)^2.*(la(i-1).*Dla(:,:,i-1)+la(1).*Dla(:,:,1))...
                                 +(8-26*la(i)).*Dla1la2la3+la1la2la3*(-26*Dla(:,:,i));
                          

                             
i = 1;
Dphip(:,:,10) = Dla1la2la3.*(24*la(i)-14)+24*la1la2la3.*Dla(:,:,i)...
                                +32*(la(3)*la(i+1)^2.*Dla(:,:,3)+la(i+1)*la(3)^2.*Dla(:,:,i+1));
                            
                            
                            
i = 2;
Dphip(:,:,11) = Dla1la2la3.*(24*la(i)-14)+24*la1la2la3.*Dla(:,:,i)...
                                +32*(la(i-1)*la(i+1)^2.*Dla(:,:,i-1)+la(i+1)*la(i-1)^2.*Dla(:,:,i+1));


i = 3;                           
Dphip(:,:,12) = Dla1la2la3.*(24*la(i)-14)+24*la1la2la3.*Dla(:,:,i)...
                                +32*(la(i-1)*la(1)^2.*Dla(:,:,i-1)+la(1)*la(i-1)^2.*Dla(:,:,1));
                            
                            
                            
                            
for i = 1:3
    Dphip(:,:,12+i)  = 32*Dla1la2la3.*(4*la(i)-1)+128*la1la2la3.*Dla(:,:,i);
end
                            
                            
                            
                            
                            
                            
i = 1;
idx = 4;


Dphip(:,:,idx) = 3*la(i)^2*(suma(:,1)-a(:,1,i)).*Dla(:,:,i)+la(i)^3.*Dsuma1 ...
                                 +2*la(i)*(a(:,1,i).*(la(3)^2+la(i+1)^2)-(a(:,1,3)*la(3)^2 ...
                                  + a(:,1,i+1)*la(i+1)^2)).*Dla(:,:,i)...
                                  +2*la(i)^2.*(a(:,1,i).*(la(3).*Dla(:,:,3)+la(i+1).*Dla(:,:,i+1))...
                                  - (a(:,1,3).*la(3).*Dla(:,:,3)+a(:,1,i+1).*la(i+1).*Dla(:,:,i+1)))...
                                  +0.25*Dla1la2la3.*(a(:,1,i).*(20*la(i)-6)+a(:,1,3).*(4*la(3)...
                                  -8*la(i)+1)+a(:,1,i+1).*(4*la(i+1)-8*la(i)+1))...
                                  +0.25*la1la2la3.*(20*a(:,1,i).*Dla(:,:,i)...
                                  +a(:,1,3).*(4*Dla(:,:,3)-8*Dla(:,:,i))...
                                  +a(:,1,i+1).*(4*Dla(:,:,i+1)-8*Dla(:,:,i)));
                              
Dphip(:,:,idx+3) = 3*la(i)^2*(suma(:,2)-a(:,2,i)).*Dla(:,:,i)+la(i)^3.*Dsuma2 ...
                                 +2*la(i)*(a(:,2,i).*(la(3)^2+la(i+1)^2)-(a(:,2,3)*la(3)^2 ...
                                  + a(:,2,i+1)*la(i+1)^2)).*Dla(:,:,i)...
                                  +2*la(i)^2.*(a(:,2,i).*(la(3).*Dla(:,:,3)+la(i+1).*Dla(:,:,i+1))...
                                  - (a(:,2,3).*la(3).*Dla(:,:,3)+a(:,2,i+1).*la(i+1).*Dla(:,:,i+1)))...
                                  +0.25*Dla1la2la3.*(a(:,2,i).*(20*la(i)-6)+a(:,2,3).*(4*la(3)...
                                  -8*la(i)+1)+a(:,2,i+1).*(4*la(i+1)-8*la(i)+1))...
                                  +0.25*la1la2la3.*(20*a(:,2,i).*Dla(:,:,i)...
                                  +a(:,2,3).*(4*Dla(:,:,3)-8*Dla(:,:,i))...
                                  +a(:,2,i+1).*(4*Dla(:,:,i+1)-8*Dla(:,:,i)));
                              
i = 2;
idx = 5;


Dphip(:,:,idx) = 3*la(i)^2*(suma(:,1)-a(:,1,i)).*Dla(:,:,i)+la(i)^3.*Dsuma1 ...
                                 +2*la(i)*(a(:,1,i).*(la(i-1)^2+la(i+1)^2)-(a(:,1,i-1)*la(i-1)^2 ...
                                  + a(:,1,i+1)*la(i+1)^2)).*Dla(:,:,i)...
                                  +2*la(i)^2.*(a(:,1,i).*(la(i-1).*Dla(:,:,i-1)+la(i+1).*Dla(:,:,i+1))...
                                  - (a(:,1,i-1).*la(i-1).*Dla(:,:,i-1)+a(:,1,i+1).*la(i+1).*Dla(:,:,i+1)))...
                                  +0.25*Dla1la2la3.*(a(:,1,i).*(20*la(i)-6)+a(:,1,i-1).*(4*la(i-1)...
                                  -8*la(i)+1)+a(:,1,i+1).*(4*la(i+1)-8*la(i)+1))...
                                  +0.25*la1la2la3.*(20*a(:,1,i).*Dla(:,:,i)...
                                  +a(:,1,i-1).*(4*Dla(:,:,i-1)-8*Dla(:,:,i))...
                                  +a(:,1,i+1).*(4*Dla(:,:,i+1)-8*Dla(:,:,i)));
                              
Dphip(:,:,idx+3) = 3*la(i)^2*(suma(:,2)-a(:,2,i)).*Dla(:,:,i)+la(i)^3.*Dsuma2 ...
                                 +2*la(i)*(a(:,2,i).*(la(i-1)^2+la(i+1)^2)-(a(:,2,i-1)*la(i-1)^2 ...
                                  + a(:,2,i+1)*la(i+1)^2)).*Dla(:,:,i)...
                                  +2*la(i)^2.*(a(:,2,i).*(la(i-1).*Dla(:,:,i-1)+la(i+1).*Dla(:,:,i+1))...
                                  - (a(:,2,i-1).*la(i-1).*Dla(:,:,i-1)+a(:,2,i+1).*la(i+1).*Dla(:,:,i+1)))...
                                  +0.25*Dla1la2la3.*(a(:,2,i).*(20*la(i)-6)+a(:,2,i-1).*(4*la(i-1)...
                                  -8*la(i)+1)+a(:,2,i+1).*(4*la(i+1)-8*la(i)+1))...
                                  +0.25*la1la2la3.*(20*a(:,2,i).*Dla(:,:,i)...
                                  +a(:,2,i-1).*(4*Dla(:,:,i-1)-8*Dla(:,:,i))...
                                  +a(:,2,i+1).*(4*Dla(:,:,i+1)-8*Dla(:,:,i)));
                              
                              
                              
                              
i = 3;
idx = 6;


Dphip(:,:,idx) = 3*la(i)^2*(suma(:,1)-a(:,1,i)).*Dla(:,:,i)+la(i)^3.*Dsuma1 ...
                                 +2*la(i)*(a(:,1,i).*(la(i-1)^2+la(1)^2)-(a(:,1,i-1)*la(i-1)^2 ...
                                  + a(:,1,1)*la(1)^2)).*Dla(:,:,i)...
                                  +2*la(i)^2.*(a(:,1,i).*(la(i-1).*Dla(:,:,i-1)+la(1).*Dla(:,:,1))...
                                  - (a(:,1,i-1).*la(i-1).*Dla(:,:,i-1)+a(:,1,1).*la(1).*Dla(:,:,1)))...
                                  +0.25*Dla1la2la3.*(a(:,1,i).*(20*la(i)-6)+a(:,1,i-1).*(4*la(i-1)...
                                  -8*la(i)+1)+a(:,1,1).*(4*la(1)-8*la(i)+1))...
                                  +0.25*la1la2la3.*(20*a(:,1,i).*Dla(:,:,i)...
                                  +a(:,1,i-1).*(4*Dla(:,:,i-1)-8*Dla(:,:,i))...
                                  +a(:,1,1).*(4*Dla(:,:,1)-8*Dla(:,:,i)));
                              
Dphip(:,:,idx+3) = 3*la(i)^2*(suma(:,2)-a(:,2,i)).*Dla(:,:,i)+la(i)^3.*Dsuma2 ...
                                 +2*la(i)*(a(:,2,i).*(la(i-1)^2+la(1)^2)-(a(:,2,i-1)*la(i-1)^2 ...
                                  + a(:,2,1)*la(1)^2)).*Dla(:,:,i)...
                                  +2*la(i)^2.*(a(:,2,i).*(la(i-1).*Dla(:,:,i-1)+la(1).*Dla(:,:,1))...
                                  - (a(:,2,i-1).*la(i-1).*Dla(:,:,i-1)+a(:,2,1).*la(1).*Dla(:,:,1)))...
                                  +0.25*Dla1la2la3.*(a(:,2,i).*(20*la(i)-6)+a(:,2,i-1).*(4*la(i-1)...
                                  -8*la(i)+1)+a(:,2,1).*(4*la(1)-8*la(i)+1))...
                                  +0.25*la1la2la3.*(20*a(:,2,i).*Dla(:,:,i)...
                                  +a(:,2,i-1).*(4*Dla(:,:,i-1)-8*Dla(:,:,i))...
                                  +a(:,2,1).*(4*Dla(:,:,1)-8*Dla(:,:,i)));
                              
                              
                              
                              
                              
                              
                              
end





end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%HermiteP3 L2���
function err = getL2H3errors(node,elemdof,elem,area,exactu,uh)
%%�õ�Hermite P3Ԫ��L2���
option.fquadorder = 6;
err = 0;
[lambda,w] = quadpts(option.fquadorder);
nQuad = size(lambda,1);
for p = 1:nQuad
        phi = getphis(lambda(p,:),node,elem);
       uhp = 0;

    for i = 1:10
        uhp = uhp +uh(elemdof(:,i)).*phi(:,i);
    end

    pxy = node(elem(:,1),:)*lambda(p,1)+node(elem(:,2),:)*lambda(p,2)+...
                       node(elem(:,3),:)*lambda(p,3);

    up = exactu(pxy);

    err = err + w(p)*(up-uhp).^2;
end

err = err.*area;
err = sqrt(sum(err));






function phi = getphis(la,node,elem)
a = zeros(size(elem,1),2,3);
a(:,:,1) = node(elem(:,1),:);
a(:,:,2) = node(elem(:,2),:);
a(:,:,3) = node(elem(:,3),:);



phi = zeros(size(elem,1),10);
la1la2la3 = la(1).*la(2).*la(3);


i = 1;
phi(:,1) = la(i).^3+3*la(i).^2.*(la(i+1)+la(3))-7*la1la2la3;
          
                     
i = 2;
phi(:,2) = la(i).^3+3*la(i).^2.*(la(i+1)+la(i-1))-7*la1la2la3;
                     
i = 3;
phi(:,3) = la(i).^3+3*la(i).^2.*(la(1)+la(i-1))-7*la1la2la3;


i = 1;
idx = 4;

phi(:,[idx,idx+3]) = la(i).^2.*(la(i+1).*(a(:,:,i+1)-a(:,:,i))+la(3).*(a(:,:,3)-a(:,:,i)))...
                                            -(a(:,:,i+1)+a(:,:,3)-2*a(:,:,i))*la1la2la3;
                                        
                                        
                                        
i = 2;
idx = 5;
phi(:,[idx,idx+3]) = la(i).^2.*(la(i+1).*(a(:,:,i+1)-a(:,:,i))+la(i-1).*(a(:,:,i-1)-a(:,:,i)))...
                                            -(a(:,:,i+1)+a(:,:,i-1)-2*a(:,:,i))*la1la2la3;
                                        
                                        
i = 3;
idx = 6;
phi(:,[idx,idx+3]) = la(i).^2.*(la(1).*(a(:,:,1)-a(:,:,i))+la(i-1).*(a(:,:,i-1)-a(:,:,i)))...
                                            -(a(:,:,1)+a(:,:,i-1)-2*a(:,:,i))*la1la2la3;
                                        
                                        
phi(:,10) = 27*la1la2la3;


% % % phi(:,idx) = 2.*area.*(la(i).^2.*(la(i+1).*e3(:,1)-la(3).*e2(:,1))...
% % %                                        +(e2(:,1)-e3(:,1)).*la1la2la3)./e_T23;
% % % phi(:,idx+3) = 2.*area.*(la(i).^2.*(la(i+1).*e3(:,2)-la(3).*e2(:,2))...
% % %                                        +(e2(:,2)-e3(:,2)).*la1la2la3)./e_T23;






end




end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%HermiteP4 L2���
function err = getL2H4errors(node,elemdof,elem,area,exactu,uh)
%%�õ�Hermite P4Ԫ��L2���
option.fquadorder = 6;
err = 0;
[lambda,w] = quadpts(option.fquadorder);
nQuad = size(lambda,1);
for p = 1:nQuad
        phi = getphis(lambda(p,:),node,elem);
       uhp = 0;

    for i = 1:15
        uhp = uhp +uh(elemdof(:,i)).*phi(:,i);
    end

    pxy = node(elem(:,1),:)*lambda(p,1)+node(elem(:,2),:)*lambda(p,2)+...
                       node(elem(:,3),:)*lambda(p,3);

    up = exactu(pxy);

    err = err + w(p)*(up-uhp).^2;
end

err = err.*area;
err = sqrt(sum(err));








function phi = getphis(la,node,elem)
a = zeros(size(elem,1),2,3);
a(:,:,1) = node(elem(:,1),:);
a(:,:,2) = node(elem(:,2),:);
a(:,:,3) = node(elem(:,3),:);



phi = zeros(size(elem,1),15);
la1la2la3 = la(1).*la(2).*la(3);
suma = la(1).*a(:,:,1)+la(2).*a(:,:,2)+la(3).*a(:,:,3);

i = 1;
phi(:,1) = 4*la(i).^3-3*la(i).^4-5*la(i).^2.*(la(i+1).^2+la(3).^2)+...
                     la1la2la3.*(8-26*la(i));
                 
          
                     
i = 2;
phi(:,2) = 4*la(i).^3-3*la(i).^4-5*la(i).^2.*(la(i+1).^2+la(i-1).^2)+...
                     la1la2la3.*(8-26*la(i));
                     
i = 3;
phi(:,3) = 4*la(i).^3-3*la(i).^4-5*la(i).^2.*(la(1).^2+la(i-1).^2)+...
                     la1la2la3.*(8-26*la(i));


                 
                 
                 
i = 1;
phi(:,10) = la1la2la3.*(24*la(i)-14)+16*la(i+1).^2.*la(3).^2;


i = 2;
phi(:,11) = la1la2la3.*(24*la(i)-14)+16*la(i+1).^2.*la(i-1).^2;


i = 3;
phi(:,12) = la1la2la3.*(24*la(i)-14)+16*la(1).^2.*la(i-1).^2;

                 
                 
for i = 1:3
  phi(:,12+i) =32* la1la2la3.*(4*la(i)-1);
end               
  


i = 1;
idx = 4;

phi(:,[idx,idx+3]) = la(i).^3.*(suma-a(:,:,i))+la(i).^2.*((la(i+1).^2 ...
                                          +la(3).^2).*a(:,:,i)-(la(i+1).^2.*a(:,:,i+1)+...
                                          la(3).^2.*a(:,:,3)))...
                                          +0.25*la1la2la3.*(a(:,:,i).*(20*la(i)-6)+a(:,:,3).*...
                                          (4*la(3)-8*la(i)+1)+a(:,:,i+1).*(4*la(i+1)-8*la(i)+1));
                                      
i = 2;
idx = 5;                                      
 
phi(:,[idx,idx+3]) = la(i).^3.*(suma-a(:,:,i))+la(i).^2.*((la(i+1).^2 ...
                                          +la(i-1).^2).*a(:,:,i)-(la(i+1).^2.*a(:,:,i+1)+...
                                          la(i-1).^2.*a(:,:,i-1)))...
                                          +la1la2la3.*(a(:,:,i).*(20*la(i)-6)+a(:,:,i-1).*...
                                          (4*la(i-1)-8*la(i)+1)+a(:,:,i+1).*(4*la(i+1)-8*la(i)+1))/4;
                                      
i = 3;
idx = 6;                                      
 
phi(:,[idx,idx+3]) = la(i).^3.*(suma-a(:,:,i))+la(i).^2.*((la(1).^2 ...
                                          +la(i-1).^2).*a(:,:,i)-(la(1).^2.*a(:,:,1)+...
                                          la(i-1).^2.*a(:,:,i-1)))...
                                          +la1la2la3.*(a(:,:,i).*(20*la(i)-6)+a(:,:,i-1).*...
                                          (4*la(i-1)-8*la(i)+1)+a(:,:,1).*(4*la(1)-8*la(i)+1))/4;
                                                                              


end




end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%���ݼ�
function pde=Possion2data2
%%%2D��Possion����
%%% u=exp(x+y)*(x^3-2)*(y^8+1);
%%% Dirichlet�߽�����



pde = struct('f', @f,  'exactu',@exactu,'exactDu',@exactDu,'g_D',@exactu,'Dg_d',@exactDu,'g_N',@g_N);

%%�Ӻ���

    function s=f(p)
       x=p(:,1); y=p(:,2);
       
       s =-2*exp(x + y).*(x.^3.*y.^8 + 8*x.^3.*y.^7 + 28*x.^3.*y.^6 ...
                      + x.^3 + 3*x.^2.*y.^8 + 3*x.^2 + 3*x.*y.^8 + 3*x - 2*y.^8 - 16*y.^7 - 56*y.^6 - 2);
      
       
    end


    function s=exactu(p)
         x=p(:,1); y=p(:,2);
         
         s = exp(x+y).*(x.^3-2).*(y.^8+1);
         
    end



    function s = exactDu(p)
         x=p(:,1); y=p(:,2);
         
         s = zeros(length(x),2);
         
         s(:,1) = exp(x + y).*(y.^8 + 1).*(x.^3 + 3*x.^2 - 2);
         s(:,2) = exp(x + y).*(x.^3 - 2).*(y.^8 + 8*y.^7 + 1);
      
    end






function s = g_N(p)
        x=p(:,1); y=p(:,2);
        
        uD = exactDu(p);
%           uD = zeros(length(x),2);
%          
%          uD(:,1) = pi*cos(pi*x).*cos(pi*y);
%          uD(:,2) = -pi*sin(pi*x).*sin(pi*y);
      
        
        s = zeros( size(p,1),1);
        
        gg = (abs(x)<1e-15);
        s(gg) =-uD(gg,1);

        
        gg = (abs(1-x)<1e-15);
        s(gg) = uD(gg,1);
       
        
        
        
        gg = (abs(y)<1e-15);
        s(gg) = -uD(gg,2);
       
        
        gg = (abs(1-y)<1e-15);
        s(gg) = uD(gg,2);
      
        
   
    end
        

end
         

