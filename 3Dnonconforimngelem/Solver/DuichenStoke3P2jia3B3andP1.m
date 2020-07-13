function  [soln,eqn] = DuichenStoke3P2jia3B3andP1(node,elem,bdFlag,pde,option)
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



if ~exist('option','var'), option = []; end


t = cputime;
%% Construct Data Structure
[elem3dof2,face] = dof3P2jia(elem);
N = size(node,1);  NT = size(elem,1);  Nu = 3*size(face,1) + NT;   Np = N;

%% Compute geometric quantities and gradient of local basis
[Dlambda,volume] = gradbasis3(node,elem);

%% Assemble stiffness matrix for Laplace operator
% generate sparse pattern
ii = zeros(91*NT,1); jj = zeros(91*NT,1); 
index = 0;
for i=1:13
    for j=i:13
         ii(index+1:index+NT) = double(elem3dof2(:,i)); 
         jj(index+1:index+NT) = double(elem3dof2(:,j));  
        index = index + NT;
    end
end

% quadrature points
if ~isfield(pde,'nu'), pde.nu = []; end
if ~isfield(option,'quadorder')
    % constant viscosity
    option.quadorder = 5;        % default order
     if ~isempty(pde.nu) && isnumeric(pde.nu) % numerical viscosity
        option.quadorder = 5;    % exact for linear diffusion coefficient
    end
end
[lambda, w] = quadpts3(option.quadorder);



A = getA(lambda,Dlambda,w,pde,node,elem,volume,NT,ii,jj,Nu);

%% Assemble the matrix for divergence operator
B=divergence(Np,Nu,lambda,Dlambda,w,elem,volume,elem3dof2);

%% Assemble right hand side
[f1,f2,f3]=righthand(Nu,option,pde,NT,node,elem,elem3dof2,volume);

%% Boundary Conditions
[AD,BD,f,u,p,g,ufreeDof,pDof] = getbdStokesP2jiaP1(node,face,bdFlag,f1,f2,f3,...
                                                                                      Np,Nu,elem3dof2,A,B,pde,elem);
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
if isempty(ufreeDof), return; end
if isempty(option) || ~isfield(option,'solver')    % no option.solver
    if length(f)+length(g) <= 1e3  % Direct solver for small size systems
        option.solver = 'direct';
    else          % Multigrid-type  solver for large size systems
      % option.solver = 'asmg';
       option.solver = 'direct';
    end
end
solver = option.solver;

% solve
switch solver
    case 'none'
        info = struct('solverTime',[],'itStep',0,'err',[],'flag',3,'stopErr',[]);        
    case 'direct'
        t = cputime;
        bigA = [AD, BD'; ...
                BD, sparse(Np,Np)];
        bigF = [f; g];
        bigu = [u; p];
        bigFreeDof = [ufreeDof; 3*Nu+pDof];        
        bigu(bigFreeDof) = bigA(bigFreeDof,bigFreeDof)\bigF(bigFreeDof);
        u = bigu(1:3*Nu);
        p = bigu(3*Nu+1:end);
        residual = norm(bigF - bigA*bigu);
        info = struct('solverTime',cputime - t,'itStep',0,'err',residual,'flag',2,'stopErr',residual);     
    case 'mg'
%         option.tol = Np^(-2);        
        option.solver  = 'WCYCLE';
        [u(ufreeDof),p,info] = mgstokes(A(ufreeDof,ufreeDof),B(:,ufreeDof),f(ufreeDof),g,...
                                        u(ufreeDof),p,elem,ufreeDof,option);         
    case 'asmg'
        [u(ufreeDof),p,info] = asmgstokes(A(ufreeDof,ufreeDof),B(:,ufreeDof),f(ufreeDof),g,...
                                          u,p,node,elem,bdFlag,ufreeDof,option); 
end

%% Post-process
if length(pDof) ~= Np % p is unique up to a constant
    % impose the condition int(p)=0
    c = sum(mean(p(elem),2).*volume)/sum(volume);
    p = p - c;
end

%% Output
soln = struct('u',u,'p',p,'Nu',Nu,'Np',Np);
eqn = struct('A',AD,'B',BD,'Lap',A,'f',f,'g',g,...
             'face',face,'ufreeDof',ufreeDof,'pDof',pDof);
info.assembleTime = assembleTime;



end







    function [AD,BD,f,u,po,g,ufreeDof,pDof] = getbdStokesP2jiaP1(node,face,bdFlag,f1,f2,f3,...
                                                                                      Np,Nu,elem3dof2,A,B,pde,elem)
    %% Boundary condition of Stokes equation: P2jia-P1 elements

    %% Initial set up
    g = zeros(Np,1);
    u = zeros(3*Nu,1);    
    po = zeros(Np,1);
    ufreeDof = (1:Nu)';
    pDof = (1:Np)';
    
    if ~exist('bdFlag','var'), bdFlag = []; end
    if ~isfield(pde,'g_D'), pde.g_D = []; end
    if ~isfield(pde,'g_N'), pde.g_N = []; end
    if ~isfield(pde,'g_R'), pde.g_R = []; end

    %% Part 1: Find Dirichlet dof and modify the matrix
    % Find Dirichlet boundary dof: fixedDof and pDof
    isFixedDof = false(Nu,1);     
    if ~isempty(bdFlag)       % case: bdFlag is not empty 
        isFixedDof(elem3dof2(bdFlag(:,1) == 1,[1,2,3])) = true;
        isFixedDof(elem3dof2(bdFlag(:,2) == 1,[4,5,6])) = true;        
        isFixedDof(elem3dof2(bdFlag(:,3) == 1,[7,8,9])) = true;
        isFixedDof(elem3dof2(bdFlag(:,4) == 1,[10,11,12])) = true;
        fixedDof = find(isFixedDof);
        ufreeDof = find(~isFixedDof);        
    end
    if isempty(fixedDof) % pure Neumann boundary condition
        % pde.g_N could be empty which is homogenous Neumann boundary condition
        fixedDof = 1;
        ufreeDof = (2:Nu)';    % eliminate the kernel by enforcing u(1) = 0;
    end
  
    % Modify the matrix
    % Build Dirichlet boundary condition into the matrix AD by enforcing
    % AD(fixedDof,fixedDof)=I, AD(fixedDof,ufreeDof)=0, AD(ufreeDof,fixedDof)=0.
    % BD(:,fixedDof) = 0 and thus BD'(fixedDof,:) = 0.
    bdidx = zeros(3*Nu,1); 
    bdidx([fixedDof; Nu+fixedDof;2*Nu+fixedDof]) = 1;
    Tbd = spdiags(bdidx,0,3*Nu,3*Nu);
    T = spdiags(1-bdidx,0,3*Nu,3*Nu);
    AD = T*A*T + Tbd;
    BD = B*T;

    %% Part 2: Find boundary faces and modify the right hand side f and g
    % Find boundary edges: Neumann 

        
    % Neumann boundary condition
    if ~isempty(pde.g_N)  
        [lambda,w] = quadpts(4);
        nQuad = size(lambda,1);
        % quadratic 

        totalface(:,:,1) = elem(:,[2,3,4]);
        totalface(:,:,2) = elem(:,[3,4,1]);
        totalface(:,:,3) = elem(:,[4,1,2]);
        totalface(:,:,4) = elem(:,[1,2,3]);
        
        
       for j = 1:4
           isNeumann = (bdFlag(:,j)==2);
           Neumannidx = elem3dof2(isNeumann,:);
          
           if ~isempty(Neumannidx)
                 
                    Neumannidx = Neumannidx(:);
                    Neumann = totalface(:,:,j);
                    Neumann = Neumann(isNeumann,:);
                    areafaces = getareaface(node,Neumann);
   
               for pp = 1:nQuad
                   lam = lambda(pp,:);
                   bdphiface = getbdphiface(lam,j);

                   pxyz = lam(1)*node(Neumann(:,1),:)+lam(2)*node(Neumann(:,2),:)...
                                      +lam(3)*node(Neumann(:,3),:);
                                     
                   gp = pde.g_N(pxyz);
                   g1 = w(pp)*[areafaces.*gp(:,1).*bdphiface(:,1);areafaces.*gp(:,1).*bdphiface(:,2); ...
                                        areafaces.*gp(:,1).*bdphiface(:,3);areafaces.*gp(:,1).*bdphiface(:,4);...
                                        areafaces.*gp(:,1).*bdphiface(:,5);areafaces.*gp(:,1).*bdphiface(:,6);...
                                        areafaces.*gp(:,1).*bdphiface(:,7);areafaces.*gp(:,1).*bdphiface(:,8);...
                                        areafaces.*gp(:,1).*bdphiface(:,9);areafaces.*gp(:,1).*bdphiface(:,10);...
                                        areafaces.*gp(:,1).*bdphiface(:,11);areafaces.*gp(:,1).*bdphiface(:,12);...
                                        areafaces.*gp(:,1).*bdphiface(:,13)];
                   
                   g2 =  w(pp)*[areafaces.*gp(:,2).*bdphiface(:,1);areafaces.*gp(:,2).*bdphiface(:,2); ...
                                          areafaces.*gp(:,2).*bdphiface(:,3);areafaces.*gp(:,2).*bdphiface(:,4);...
                                          areafaces.*gp(:,2).*bdphiface(:,5);areafaces.*gp(:,2).*bdphiface(:,6);...
                                          areafaces.*gp(:,2).*bdphiface(:,7);areafaces.*gp(:,2).*bdphiface(:,8);...
                                          areafaces.*gp(:,2).*bdphiface(:,9);areafaces.*gp(:,2).*bdphiface(:,10);...
                                          areafaces.*gp(:,2).*bdphiface(:,11);areafaces.*gp(:,2).*bdphiface(:,12);...
                                          areafaces.*gp(:,2).*bdphiface(:,13)];
                   
                   g3 =  w(pp)*[areafaces.*gp(:,3).*bdphiface(:,1);areafaces.*gp(:,3).*bdphiface(:,2); ...
                                          areafaces.*gp(:,3).*bdphiface(:,3);areafaces.*gp(:,3).*bdphiface(:,4);...
                                          areafaces.*gp(:,3).*bdphiface(:,5);areafaces.*gp(:,3).*bdphiface(:,6);...
                                          areafaces.*gp(:,3).*bdphiface(:,7);areafaces.*gp(:,3).*bdphiface(:,8);...
                                          areafaces.*gp(:,3).*bdphiface(:,9);areafaces.*gp(:,3).*bdphiface(:,10);...
                                          areafaces.*gp(:,3).*bdphiface(:,11);areafaces.*gp(:,3).*bdphiface(:,12);...
                                          areafaces.*gp(:,3).*bdphiface(:,13)];
                        
                                  
                    S = size(f1,1);
                    f1 = f1 + accumarray(Neumannidx,g1,[S,1]);
                    f2 = f2 + accumarray(Neumannidx,g2,[S,1]);
                    f3 = f3 + accumarray(Neumannidx,g3,[S,1]);


                   
                   
                end
            end
        end
    end                       
if isempty(Neumannidx)
                 Neumann = [];
end
       

        
    f = [f1; f2;f3];
    % The case non-empty Neumann but g_N=[] corresponds to the zero flux
    % boundary condition on Neumann edges and no modification is needed.

    % Dirichlet boundary conditions
    if ~isempty(fixedDof) && ~isempty(pde.g_D) && ~(isnumeric(pde.g_D) && (pde.g_D == 0))
         NF = length(face);
        idx = (fixedDof<=NF);
        isDirichlet = fixedDof(idx);
        Dirichletface = face(isDirichlet,:);
        u1 = zeros(Nu,1);
        u2 = zeros(Nu,1);
        u3 = zeros(Nu,1);           
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uD = pde.g_D(node(fixedDof,:));  % bd value at vertex dofs.................................
        %找边界面上的基函数,并对其赋值,用三角形上的数值积分
        %用5点积分来算
        [lambda,w] = quadpts(5);
        nQuad=length(lambda);
      
        for p=1:nQuad
             pxyz=lambda(p,1)*node(Dirichletface(:,1),:)+lambda(p,2)*node(Dirichletface(:,2),:)...
                                +lambda(p,3)*node(Dirichletface(:,3),:);
             uD=pde.g_D(pxyz);
             uD1 = uD(:,1);
             uD2 = uD(:,2);
             uD3 = uD(:,3);
             
             u1D = w(p)*[uD1*lambda(p,1);uD1*lambda(p,2);uD1*lambda(p,3)];
             u2D = w(p)*[uD2*lambda(p,1);uD2*lambda(p,2);uD2*lambda(p,3)];
             u3D = w(p)*[uD3*lambda(p,1);uD3*lambda(p,2);uD3*lambda(p,3)];

             u1(fixedDof)=u1(fixedDof)+u1D;
             u2(fixedDof)=u2(fixedDof)+u2D;
             u3(fixedDof)=u3(fixedDof)+u3D;
             
        end
        
        
        u = [u1; u2;u3]; % Dirichlet bd condition is built into u
        f = f - A*u;  % bring affect of nonhomgenous Dirichlet bd condition to
        g = g - B*u;  % the right hand side
        g = g - mean(g);         
        f(fixedDof) = u1(fixedDof);
        f(fixedDof+Nu) = u2(fixedDof);
        f(fixedDof+2*Nu) = u3(fixedDof);
    end
    % The case non-empty Dirichlet but g_D=[] corresponds to the zero Dirichlet
    % boundary condition and no modification is needed.
    
    % modfiy pressure dof for pure Dirichlet
    if isempty(Neumann)
        pDof = (1:Np-1)';
    end
    

    
    ufreeDof = [ufreeDof; Nu+ufreeDof;2*Nu+ufreeDof];
    end 




function A =  getA(lambda,Dlambda,w,pde,node,elem,volume,NT,ii,jj,Nu)
nQuad = size(lambda,1);

sA11=zeros(91*NT,nQuad);
sA12=zeros(91*NT,nQuad);
sA13=zeros(91*NT,nQuad);
sA21=zeros(91*NT,nQuad);
sA22=zeros(91*NT,nQuad);
sA23=zeros(91*NT,nQuad);
sA31=zeros(91*NT,nQuad);
sA32=zeros(91*NT,nQuad);
sA33=zeros(91*NT,nQuad);
sA = zeros(91*NT,nQuad);

for pp=1:nQuad

    Dphip = getDphipP2jia(Dlambda,lambda(pp,:));

     index = 0;
    for i = 1:13
        for j = i:13
            aij = 0;
            a11ij = 0;
            a12ij = 0;
            a13ij = 0;
            a21ij = 0;
            a22ij = 0;
            a23ij = 0;
            a31ij = 0;
            a32ij = 0;
            a33ij = 0;
            
            Dphipi = Dphip(:,:,i);
            Dphipj = Dphip(:,:,j);
            a = zeros(NT,3,3);
            for kk = 1:3
                for ll = 1:3
                     a(:,kk,ll) = Dphipi(:,ll).*Dphipj(:,kk);
                end
            end

            
              if isempty(pde.nu) || isnumeric(pde.nu)
                %Aij = Aij + w(pp)*dot(Dphip(:,:,i),Dphip(:,:,j),2);
                aij = aij+w(pp)*(a(:,1,1)+a(:,2,2)+a(:,3,3));
                a11ij = a11ij + w(pp)*a(:,1,1);
                a12ij = a12ij + w(pp)*a(:,1,2);
                a13ij = a13ij + w(pp)*a(:,1,3);
                a21ij = a21ij + w(pp)*a(:,2,1);
                a22ij = a22ij + w(pp)*a(:,2,2);
                a23ij = a23ij + w(pp)*a(:,2,3);
                a31ij = a31ij + w(pp)*a(:,3,1);
                a32ij = a32ij + w(pp)*a(:,3,2);
                a33ij = a33ij + w(pp)*a(:,3,3);
            else
                pxyz = lambda(pp,1)*node(elem(:,1),:) ...
                    + lambda(pp,2)*node(elem(:,2),:) ...
                    + lambda(pp,3)*node(elem(:,3),:)...
                    +lambda(pp,4)*node(elem(:,4),:);
                 % Aij = Aij + w(pp)*dot(Dphip(:,:,i),Dphip(:,:,j),2).*pde.nu(pxyz);
                 nu = pde.nu(pxyz);
                aij = aij + w(pp)*(a(:,1,1)+a(:,2,2)+a(:,3,3)).*nu;
                a11ij = a11ij + w(pp)*a(:,1,1).*nu;
                a12ij = a12ij + w(pp)*a(:,1,2).*nu;
                a13ij = a13ij + w(pp)*a(:,1,3).*nu;
                a21ij = a21ij + w(pp)*a(:,2,1).*nu;
                a22ij = a22ij + w(pp)*a(:,2,2).*nu;
                a23ij = a23ij + w(pp)*a(:,2,3).*nu;
                a31ij = a31ij + w(pp)*a(:,3,1).*nu;
                a32ij = a32ij + w(pp)*a(:,3,2).*nu;
                a33ij = a33ij + w(pp)*a(:,3,3).*nu;
              end
              if ~isempty(pde.nu) && (pde.nu~=1)
                    %Aij = pde.nu*Aij;
                aij = aij + w(pp)*(a(:,1,1)+a(:,2,2)+a(:,3,3)).*pde.nu;
                a11ij = a11ij + w(pp)*a(:,1,1).*pde.nu;
                a12ij = a12ij + w(pp)*a(:,1,2).*pde.nu;
                a13ij = a13ij + w(pp)*a(:,1,3).*pde.nu;
                a21ij = a21ij + w(pp)*a(:,2,1).*pde.nu;
                a22ij = a22ij + w(pp)*a(:,2,2).*pde.nu;
                a23ij = a23ij + w(pp)*a(:,2,3).*pde.nu;
                a31ij = a31ij + w(pp)*a(:,3,1).*pde.nu;
                a32ij = a32ij + w(pp)*a(:,3,2).*pde.nu;
                a33ij = a33ij + w(pp)*a(:,3,3).*pde.nu;
             end
           % Aij = Aij.*volume;
                aij = aij.*volume;
                a11ij = a11ij.*volume;
                a12ij = a12ij.*volume;
                a13ij = a13ij.*volume;
                a21ij = a21ij.*volume;
                a22ij = a22ij.*volume;
                a23ij = a23ij.*volume;
                a31ij = a31ij.*volume;
                a32ij = a32ij.*volume;
                a33ij = a33ij.*volume;
            sA(index+1:index+NT,pp) = aij;
            sA11(index+1:index+NT,pp) = aij+a11ij;
            sA12(index+1:index+NT,pp) = a12ij;
            sA13(index+1:index+NT,pp) = a13ij;
            sA21(index+1:index+NT,pp) = a21ij;
            sA22(index+1:index+NT,pp) =aij+ a22ij;
            sA23(index+1:index+NT,pp) = a23ij;
            sA31(index+1:index+NT,pp) = a31ij;
            sA32(index+1:index+NT,pp) = a32ij;
            sA33(index+1:index+NT,pp) = aij+a33ij;
            
            index = index + NT;
        end
    end
end
     %sA = sum(sA,2);
     sA11 = sum(sA11,2);
     sA12 = sum(sA12,2);
     sA13 = sum(sA13,2);
     sA21 = sum(sA21,2);
     sA22 = sum(sA22,2);
     sA23 = sum(sA23,2);
     sA31 = sum(sA31,2);
     sA32 = sum(sA32,2);
     sA33 = sum(sA33,2);
     
     
     diagIdx = (ii == jj);   upperIdx = ~diagIdx;
     A = sparse(ii(diagIdx),jj(diagIdx),sA11(diagIdx),Nu,Nu);
     AU = sparse(ii(upperIdx),jj(upperIdx),sA11(upperIdx),Nu,Nu);
     A11 =A + AU + AU';

     A = sparse(ii(diagIdx),jj(diagIdx),sA22(diagIdx),Nu,Nu);
     AU = sparse(ii(upperIdx),jj(upperIdx),sA22(upperIdx),Nu,Nu);
     A22 =A + AU + AU';


    A = sparse(ii(diagIdx),jj(diagIdx),sA33(diagIdx),Nu,Nu);
    AU = sparse(ii(upperIdx),jj(upperIdx),sA33(upperIdx),Nu,Nu);
    A33 =A + AU + AU';


    A = sparse(ii(diagIdx),jj(diagIdx),sA12(diagIdx),Nu,Nu);
    AU = sparse(ii(upperIdx),jj(upperIdx),sA12(upperIdx),Nu,Nu);
    AD = sparse(ii(upperIdx),jj(upperIdx),sA21(upperIdx),Nu,Nu);
    A12 = A+AU+AD';


    A = sparse(ii(diagIdx),jj(diagIdx),sA13(diagIdx),Nu,Nu);
    AU = sparse(ii(upperIdx),jj(upperIdx),sA13(upperIdx),Nu,Nu);
    AD = sparse(ii(upperIdx),jj(upperIdx),sA31(upperIdx),Nu,Nu);
    A13 = A+AU+AD';


    A = sparse(ii(diagIdx),jj(diagIdx),sA23(diagIdx),Nu,Nu);
    AU = sparse(ii(upperIdx),jj(upperIdx),sA23(upperIdx),Nu,Nu);
    AD = sparse(ii(upperIdx),jj(upperIdx),sA32(upperIdx),Nu,Nu);
    A23 = A+AU+AD';
    
%     A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),Nu,Nu);
%     AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),Nu,Nu);
%     Ao = A+AU+AU';


    A = [A11,A12,A13;
        A12',A22,A23;
        A13',A23',A33];
end
 

function B=divergence(Np,Nu,lambda,Dlambda,w,elem,volume,elem3dof2)
%%计算散度项B
Dx = sparse(Np,Nu);
Dy = sparse(Np,Nu);
Dz = sparse(Np,Nu);
nQuad = size(lambda,1);



for p = 1:nQuad
    % Dphi at quadrature points


   Dphip = getDphipP2jia(Dlambda,lambda(p,:));

    
       for i = 1:13 
             for j = 1:4
                 
                   Dxij =  w(p)*Dphip(:,1,i)*lambda(p,j);
                   Dyij =  w(p)*Dphip(:,2,i)*lambda(p,j);
                   Dzij =  w(p)*Dphip(:,3,i)*lambda(p,j);
                 
                   Dx = Dx + sparse(elem(:,j),double(elem3dof2(:,i)),Dxij.*volume,Np,Nu);
                   Dy = Dy + sparse(elem(:,j),double(elem3dof2(:,i)),Dyij.*volume,Np,Nu);
                   Dz = Dz + sparse(elem(:,j),double(elem3dof2(:,i)),Dzij.*volume,Np,Nu);
                   
             end
      end
end
 B = -[Dx Dy Dz];
end



function [f1,f2,f3]=righthand(Nu,option,pde,NT,node,elem,elem3dof2,volume)
%% 计算右端项
f1 = zeros(Nu,1);
f2 = zeros(Nu,1);
f3 = zeros(Nu,1);
if ~isfield(option,'fquadorder')
    option.fquadorder = 5;   % default order
end
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if ~isempty(pde.f) 
    % quadrature points in the barycentric coordinate
    [lambda,weight] = quadpts3(option.fquadorder);
    % basis values at quadrature points
   
     phi = getphiP2jia(lambda);
     
     nQuad = size(lambda,1);
    ft1 = zeros(NT,13);
    ft2 = zeros(NT,13);
    ft3 = zeros(NT,13);
    for p = 1:nQuad
        % quadrature points in the x-y-z coordinate
        pxyz =    lambda(p,1)*node(elem(:,1),:) ...
                    + lambda(p,2)*node(elem(:,2),:) ...
                    + lambda(p,3)*node(elem(:,3),:)...
                    + lambda(p,4)*node(elem(:,4),:);
        % function values at quadrature points
        fp = pde.f(pxyz);
        % evaluate fp outside.
        for j = 1:13
            ft1(:,j) = ft1(:,j) + fp(:,1).*phi(p,j)*weight(p);
            ft2(:,j) = ft2(:,j) + fp(:,2).*phi(p,j)*weight(p);
            ft3(:,j) = ft3(:,j) + fp(:,3).*phi(p,j)*weight(p);
        end
    end
    ft1 = ft1.*repmat(volume,1,13);
    ft2 = ft2.*repmat(volume,1,13);
    ft3 = ft3.*repmat(volume,1,13);
    f1 = accumarray(elem3dof2(:),ft1(:),[Nu 1]);
    f2 = accumarray(elem3dof2(:),ft2(:),[Nu 1]);
    f3 = accumarray(elem3dof2(:),ft3(:),[Nu 1]);
end
end


   


function bdphiface = getbdphiface(lambda,j)

 n = size(lambda,1);
 la = zeros(n,1);
switch j
    case 1
        la2 = lambda(:,1);
        la3 = lambda(:,2);
        la4 = lambda(:,3);
        lam = [la,la2,la3,la4];
   
    case 2
        la3 = lambda(:,1);
        la4 = lambda(:,2);
        la1 = lambda(:,3);
        lam = [la1,la,la3,la4];
     
    case 3
        la4 = lambda(:,1);
        la1 = lambda(:,2);
        la2 = lambda(:,3);
        
        lam = [la1,la2,la,la4];
    
        
        
        
    case 4
        la1 = lambda(:,1);
        la2 = lambda(:,2);
        la3 = lambda(:,3);
        
        lam = [la1,la2,la3,la];
    
end
bdphiface = getphiP2jia(lam);

end




function areaface = getareaface(node,Neumann)
%%计算面的面积


    v12 = node(Neumann(:,2),:)-node(Neumann(:,1),:);
    v13 = node(Neumann(:,3),:)-node(Neumann(:,1),:);
    
    
    areaface = 0.5*sqrt(abs(sum(mycross(v12,v13,2).^2,2)));

end

