function  [soln,eqn] = Poisson2DP3N(node,elem,bdFlag,pde,option)
%%2D上P3非协调元求解Poisson方程
% 
%       -div(grad u)  = f in \Omega,
%   with 
%       Dirichlet boundary condition        u = g_D  on \Gamma_D, 
%       Neumann boundary condition du/dn  = g_N  on \Gamma_N.



if ~exist('option','var'), option = []; end


t = cputime;
%% Construct Data Structure
[elem2Ddof3,edge] = dof2DP3N(elem);
%N = size(node,1);
NT = size(elem,1);  
NF = size(edge,1);
Nu = 3*NF + NT;
%% Compute geometric quantities and gradient of local basis
[Dlambda,volume] = gradbasis(node,elem);

%% Assemble stiffness matrix for Laplace operator
% generate sparse pattern
ii = zeros(55*NT,1); jj = zeros(55*NT,1); 
index = 0;
for i=1:10
    for j=i:10
         ii(index+1:index+NT) = double(elem2Ddof3(:,i)); 
         jj(index+1:index+NT) = double(elem2Ddof3(:,j));  
        index = index + NT;
    end
end

% quadrature points
if ~isfield(pde,'nu'), pde.nu = []; end
if ~isfield(option,'quadorder')
    % constant viscosity
    option.quadorder = 6;        % default order
     if ~isempty(pde.nu) && isnumeric(pde.nu) % numerical viscosity
        option.quadorder = 6;    % exact for linear diffusion coefficient
    end
end

[lambda, w] = quadpts(option.quadorder);

% compute non-zeros
sA=Dphips(lambda,Dlambda,w,pde,node,elem,volume,NT);

% assemble the matrix
diagIdx = (ii == jj);   upperIdx = ~diagIdx;
A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),Nu,Nu);
AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),Nu,Nu);
A = A + AU + AU';
clear Aij ii jj sA


%% Assemble right hand side
b=righthand(Nu,option,pde,NT,node,elem,elem2Ddof3,volume);

%% Boundary Conditions
[AD,u,ufreeDof,b,isPureNeumann] = getbdPossion2DP3N(node,edge,bdFlag,b,...
                                                                                      Nu,elem2Ddof3,A,pde,elem);
%% Record assembeling time
assembleTime = cputime - t;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end


%% Record assembeling time
assembleTime = cputime - t;
if ~isfield(option,'printlevel'), option.printlevel = 1; end
if option.printlevel >= 2
    fprintf('Time to assemble matrix equation %4.2g s\n',assembleTime);
end

%% Solve the system of linear equations


u(ufreeDof) = AD(ufreeDof,ufreeDof)\b(ufreeDof);

if isPureNeumann
    intguh = double(sparse(NT,1)); 
    [lambda,w] = quadpts(3);  %This is P3 element.
    nQuad = size(lambda,1);
    phi = getphi2DP3N(lambda);
    UD = u(elem2Ddof3);
    for  p = 1:nQuad
          intguh(:) = intguh(:) + sum(w(p)*phi(p,:).*UD,2);
    end
    intguh = intguh.*volume;
    uc = sum(intguh)/sum(volume);
    u = u - uc;    % int u = 0
    
end
%% Output
soln = struct('u',u,'Nu',Nu);
eqn = struct('A',AD,'Lap',A,'f',b,...
             'edge',edge,'ufreeDof',ufreeDof);


end




    function [AD,u,ufreeDof,b,isPureNeumann] = getbdPossion2DP3N(node,edge,bdFlag,b,...
                                                                                      Nu,elem2Ddof3,A,pde,elem)
    %% Boundary condition of Stokes equation: P2jia-P1 elements

    %% Initial set up
    % set in Neumann boundary condition
    u = zeros(Nu,1);    
    if ~isfield(pde,'g_D'), pde.g_D = []; end
    if ~isfield(pde,'g_N'), pde.g_N = []; end
  
    %% Part 1: Find Dirichlet dof and modify the matrix
    % Find Dirichlet boundary dof: fixedDof and pDof
    isFixedDof = false(Nu,1);     
    if ~isempty(bdFlag)       % case: bdFlag is not empty 
        isFixedDof(elem2Ddof3(bdFlag(:,1) == 1,1:3)) = true;
        isFixedDof(elem2Ddof3(bdFlag(:,2) == 1,4:6)) = true;        
        isFixedDof(elem2Ddof3(bdFlag(:,3) == 1,7:9)) = true;
        fixedDof = find(isFixedDof);
        ufreeDof = find(~isFixedDof);        
    end
    isPureNeumann = false;
    if isempty(fixedDof) % pure Neumann boundary condition
        % pde.g_N could be empty which is homogenous Neumann boundary condition
        isPureNeumann = true;
        fixedDof = 1;
        ufreeDof = 2:Nu;    % eliminate the kernel by enforcing u(1) = 0;
    end
  
    % Modify the matrix
    % Build Dirichlet boundary condition into the matrix AD by enforcing
    % AD(fixedDof,fixedDof)=I, AD(fixedDof,ufreeDof)=0, AD(ufreeDof,fixedDof)=0.
    % BD(:,fixedDof) = 0 and thus BD'(fixedDof,:) = 0.
    Ndof = Nu;
    if ~isempty(fixedDof)
        bdidx = zeros(Ndof,1); 
        bdidx(fixedDof) = 1;
        Tbd = spdiags(bdidx,0,Ndof,Ndof);
        T = spdiags(1-bdidx,0,Ndof,Ndof);
        AD = T*A*T + Tbd;
    else
        AD = A;
    end
 
    %% Part 2: Find boundary faces and modify the right hand b
    % Find boundary faces bdFace for Neumann boundary condition
    % Neumann boundary condition
    if ~isempty(pde.g_N)  
        [lambda,w] = quadpts1(6);
        nQuad = size(lambda,1); 

        totaledge(:,:,1) = elem(:,[2,3]);
        totaledge(:,:,2) = elem(:,[3,1]);
        totaledge(:,:,3) = elem(:,[1,2]);
        
        
       for j = 1:3
           isNeumann = (bdFlag(:,j)==2);
           Neumannidx = elem2Ddof3(isNeumann,:);
          
           if ~isempty(Neumannidx)
                   Neumannidx = Neumannidx(:);
                    Neumann = totaledge(:,:,j);
                    Neumann = Neumann(isNeumann,:);
                    areaedge = getareaedge(node,Neumann);
                    g = 0;
   
               for pp = 1:nQuad
                   lam = lambda(pp,:);
                   bdphiedge = getbdphiface(lam,j);

                   pxyz = lam(1)*node(Neumann(:,1),:)+lam(2)*node(Neumann(:,2),:);
                                     
                   gp = pde.g_N(pxyz);
                        g = g+ w(pp)*areaedge.*gp.*bdphiedge;

    
               end
                g = g(:);
                S = size(b,1);
                b = b + accumarray(Neumannidx,g,[S,1]);
                
            end
        end
    end                       
                                                                   
       

        
    % The case non-empty Neumann but g_N=[] corresponds to the zero flux
    % boundary condition on Neumann edges and no modification is needed.

    % Dirichlet boundary conditions
    if ~isPureNeumann && ~isempty(fixedDof) && ...
       ~isempty(pde.g_D) && ~(isnumeric(pde.g_D) && all(pde.g_D == 0))    % nonzero g_D
                NF = length(edge);
                idx = (fixedDof<=NF);
                isDirichlet = fixedDof(idx);
                Dirichletedge = edge(isDirichlet,:);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uD = pde.g_D(node(fixedDof,:));  % bd value at vertex dofs.................................
        %找边界边上的基函数,并对其赋值,用线段上的数值积分
        [lambda,w] = quadpts1(6);
        nQuad=length(lambda);
      
        for p=1:nQuad
            pxyz=lambda(p,1)*node(Dirichletedge(:,1),:)+lambda(p,2)*node(Dirichletedge(:,2),:);
             uD=pde.g_D(pxyz);
             
            la1 = lambda(p,1);
            la2 = lambda(p,2);
             

            uD = w(p)*[uD*la1.^2;uD*la2.^2;uD*la1*la2];
            u(fixedDof)=u(fixedDof)+uD;
             
        end
        
        
        b = b - A*u;  % bring affect of nonhomgenous Dirichlet bd condition to       
    end
    % The case non-empty Dirichlet but g_D=[] corresponds to the zero Dirichlet
    % boundary condition and no modification is needed.
    
    if ~isPureNeumann % non-empty Dirichlet boundary condition
        b(fixedDof) = u(fixedDof);
    end
    
    
   if isPureNeumann
       e1 = elem2Ddof3(:,[1,2,4,5,7,8]);
       e1 = e1(:);
       e2 = elem2Ddof3(:,[3,6,9]);
       e2 = e2(:);
       
        s = sum(b(e1))/3 + sum(b(e2))/6 +sum(b(elem2Ddof3(:,10)));
        s = s/Nu;
        b = b - s;   % compatilbe condition: sum(b) = 0
        b(1) = 0;          % 1 is fixedNode and set u(1) = 0
    end

    end 



function   sA=Dphips(lambda,Dlambda,w,pde,node,elem,volume,NT)
%%%计算矩阵sA
nQuad=size(lambda,1);
sA=zeros(55*NT,nQuad);
Da1=Dlambda(:,:,1);
Da2=Dlambda(:,:,2);
Da3=Dlambda(:,:,3);

for p=1:nQuad
    la1 = lambda(p,1);
    la2 = lambda(p,2);
    la3 = lambda(p,3);
    Dphip = getDphip2DP3N(Da1,Da2,Da3,la1,la2,la3);

     index = 0;
    for i = 1:10
        for j = i:10
            Aij = 0;
              if isempty(pde.nu) || isnumeric(pde.nu)
                Aij = Aij + w(p)*dot(Dphip(:,:,i),Dphip(:,:,j),2);
            else
                pxy = lambda(p,1)*node(elem(:,1),:) ...
                    + lambda(p,2)*node(elem(:,2),:) ...
                    + lambda(p,3)*node(elem(:,3),:);
                  Aij = Aij + w(p)*dot(Dphip(:,:,i),Dphip(:,:,j),2).*pde.nu(pxy);
              end
              if ~isempty(pde.nu) && (pde.nu~=1)
                    Aij = pde.nu*Aij;
             end
            Aij = Aij.*volume;
            sA(index+1:index+NT,p) = Aij;
            index = index + NT;
        end
    end
end
     sA = sum(sA,2);
end

 




function b=righthand(Nu,option,pde,NT,node,elem,elem2Ddof3,volume)
%% 计算右端项 b
if ~isfield(option,'fquadorder')
    option.fquadorder = 5;   % default order
end


    % quadrature points in the barycentric coordinate
    [lambda,weight] = quadpts(option.fquadorder);
    % basis values at quadrature points
   
     phi = getphi2DP3N(lambda);
     
     nQuad = size(lambda,1);
    ft = zeros(NT,10);
    for p = 1:nQuad
        % quadrature points in the x-y-z coordinate
        pxyz =    lambda(p,1)*node(elem(:,1),:) ...
                    + lambda(p,2)*node(elem(:,2),:) ...
                    + lambda(p,3)*node(elem(:,3),:);
        % function values at quadrature points
        fp = pde.f(pxyz);
        % evaluate fp outside.
        for j = 1:10
            ft(:,j) = ft(:,j) + fp.*phi(p,j)*weight(p);
        end
    end
    ft = ft.*repmat(volume,1,10);
    b = accumarray(elem2Ddof3(:),ft(:),[Nu 1]);
end





function bdphiedge = getbdphiface(lambda,j)

 n = size(lambda,1);
 la = zeros(n,1);
switch j
    case 1
        la2 = lambda(:,1);
        la3 = lambda(:,2);
        lam = [la,la2,la3];
   
    case 2
        la3 = lambda(:,1);
        la1 = lambda(:,2);
        lam = [la1,la,la3];
     
    case 3
        la1 = lambda(:,1);
        la2 = lambda(:,2);      
        lam = [la1,la2,la];
    
       
    
end
bdphiedge = getphi2DP3N(lam);

end




function areaedge = getareaedge(node,Neumann)
%%计算面的面积


    v12 = node(Neumann(:,2),:)-node(Neumann(:,1),:);
    
    
    areaedge = sqrt(sum(v12.^2,2));

end


