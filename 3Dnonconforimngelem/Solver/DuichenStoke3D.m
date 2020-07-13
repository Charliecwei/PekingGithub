function [soln,eqn] = DuichenStoke3D(node,elem,bdFlag,pde,option)
%% Solve the symmetry stokes equation
% 
%       -div[2*mu*eplison(u)] + grad p = f in \Omega,
%                        - div u = 0  in \Omega,
%   with 
%       Dirichlet boundary condition        u = g_D  on \Gamma_D, 
%       Neumann boundary condition sigma  = g_N  on \Gamma_N. 
% 
%     eplison(u) = 0.5*[grad(u)+grad(u)^T]
%     sigma = (2*mu*eplison(u)- p*I)*n
%
% option deter determine the solver type

%% Construct  Data A and B
%  construct elemdof of vecter u
switch option.usolver
    case 'P2jia'
       [elemdof,face] = dof3P2jia(elem);
        quadorder = 5;
       eleminforming = struct('face',face);
        
    case '3DP2N'
       [elemdof,edge] = dof3DP2N(elem);
       quadorder = 3;
       eleminforming = struct('edge',edge);
       
    case '3DP3N'
        [elemdof,face,edge] = dof3DP3N(elem);
        quadorder = 5;
        eleminforming = struct('face',face,'edge',edge);
        
    case 'P3jia5B4'
        [elemdof,face] = dof3P3jia5B4s(elem);
        quadorder = 6;
         eleminforming = struct('face',face);
         
    case '3DP2'
        [elemdof,edge] = dof3P2(elem);
         quadorder = 3;
         eleminforming = struct('edge',edge);
         
    case '3DP3'
        [elemdof,face,edge] = dof3P3(elem);
        quadorder = 5;
        eleminforming = struct('face',face,'edge',edge);
    case '3DP3jia8B4'
         [elemdof,face] = dof3P3jia8B4s(elem);
        quadorder = 6;
         eleminforming = struct('face',face);
     
end

switch option.psolver
    case 'P1'
        elemdofp = elem;
    case 'P2'
        [elemdofp,~] = dof3P2(elem); 
end

% Compute geometric quantities and gradient of local basis
[Dlambda,volume] = gradbasis3(node,elem);

NT = size(elem,1); % The number of elem;
elemsize = size(elemdof,2); %The number of dof on a elem
Nu = double(max(elemdof(:))); % The number of u_1
Np = double(max(elemdofp(:))); % The dimension of pressure Q

% Assemble stiffness matrix for Laplace operator
% generate sparse pattern
Ngs = (elemsize+1)*elemsize/2;
ii = zeros(Ngs*NT,1); jj = zeros(Ngs*NT,1);
index = 0;        

for i=1:elemsize
    for j=i:elemsize
         ii(index+1:index+NT) = double(elemdof(:,i)); 
         jj(index+1:index+NT) = double(elemdof(:,j));  
        index = index + NT;
    end
end


if ~isfield(pde,'nu'), pde.nu = []; end
% quadrature points
[lambda, w] = quadpts3(quadorder);

%get the matrix A
A = getA(lambda,Dlambda,w,pde,node,elem,volume,NT,ii,jj,Nu,Ngs,elemsize,option);

%get the matrix B
B=divergence(Np,Nu,lambda,Dlambda,w,elemdofp,volume,elemdof,elemsize,option);

%(Np,Nu,lambda,Dlambda,w,elem3dofp2,volume,elem3dof3)
%%

%% Assemble right hand side
[f1,f2,f3]=righthand(Nu,option,pde,NT,node,elem,elemdof,volume,elemsize);

%% Boundary Conditions
[AD,BD,f,u,p,g,ufreeDof,pDof] = getbdStokes(node,eleminforming,bdFlag,f1,f2,f3,...
                                                                                      Np,Nu,elemdof,A,B,pde,elem,option);
                                                                                  
                                                                                  
                                                                                  
                                                                                  
                                                                                  
        bigA = [AD, BD'; ...
                 BD, sparse(Np,Np)];
        bigF = [f; g];
        bigu = [u; p];
        bigFreeDof = [ufreeDof; 3*Nu+pDof];        
        bigu(bigFreeDof) = bigA(bigFreeDof,bigFreeDof)\bigF(bigFreeDof);
        u = bigu(1:3*Nu);
        p = bigu(3*Nu+1:end);
     %   residual = norm(bigF - bigA*bigu);
        
    %% Post-process       
    if length(pDof) ~= Np % p is unique up to a constant
        % impose the condition int(p)=0
        intguh =0; 
        switch option.psolver
            case 'P1'
                [lambda,w] = quadpts3(1); %This is P1 element.
                 phipp = lambda;
            case 'P2'
                [lambda,w] = quadpts3(2);  %This is P2 element.
                 phipp = getphi3DP2(lambda);
        end
        nQuad = size(lambda,1); 
        UD = u(elemdofp);
        for  pp = 1:nQuad
              intguh = intguh + sum(w(pp)*phipp(pp,:).*UD,2);
        end
        intguh = intguh.*volume;
        c = sum(intguh)/sum(volume);
        p = p - c;
    end
    
    soln = struct('u',u,'p',p,'Nu',Nu,'Np',Np);
    eqn = struct('A',AD,'B',BD,'Lap',A,'f',f,'g',g,...
                    'informing',eleminforming,'ufreeDof',ufreeDof,'pDof',pDof);
end




function A =  getA(lambda,Dlambda,w,pde,node,elem,volume,NT,ii,jj,Nu,Ngs,elemsize,option)
nQuad = size(lambda,1);
sA11=zeros(Ngs*NT,nQuad);
sA12=zeros(Ngs*NT,nQuad);
sA13=zeros(Ngs*NT,nQuad);
sA21=zeros(Ngs*NT,nQuad);
sA22=zeros(Ngs*NT,nQuad);
sA23=zeros(Ngs*NT,nQuad);
sA31=zeros(Ngs*NT,nQuad);
sA32=zeros(Ngs*NT,nQuad);
sA33=zeros(Ngs*NT,nQuad);
sA = zeros(Ngs*NT,nQuad);

for pp=1:nQuad

    switch option.usolver
        case 'P2jia'
            Dphip = getDphipP2jia(Dlambda,lambda(pp,:));
            
        case '3DP2N'
              Dphip = getDphip3DP2N(Dlambda,lambda(pp,:));
              
        case '3DP3N'
               Dphip = getDphip3DP3N(Dlambda,lambda(pp,:));
               
        case 'P3jia5B4'
             Dphip = getDphipP3jiapos(Dlambda,lambda(pp,:));
             
        case '3DP2'
            Dphip = getDphip3DP2(Dlambda,lambda(pp,:));
            
        case '3DP3'
            Dphip = getDphip3DP3(Dlambda,lambda(pp,:));
            
        case '3DP3jia8B4'
             Dphip = getDphipP3jia8B4(Dlambda,lambda(pp,:));
            

            
               
    end

     index = 0;
    for i = 1:elemsize
        for j = i:elemsize
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
 


function B=divergence(Np,Nu,lambda,Dlambda,w,elemdofp,volume,elemdof,elemsize,option)

switch option.psolver
    case 'P1'
        phipp = lambda;
    case 'P2'
        phipp = getphi3DP2(lambda);
end
elemsizep = size(elemdofp,2);

Dx = sparse(Np,Nu);
Dy = sparse(Np,Nu);
Dz = sparse(Np,Nu);
nQuad = size(lambda,1);

for p = 1:nQuad
    % Dphi at quadrature points

   switch option.usolver
       case 'P2jia'
            Dphip = getDphipP2jia(Dlambda,lambda(p,:));
       case '3DP2N'
           Dphip = getDphip3DP2N(Dlambda,lambda(p,:));
       case '3DP3N'
            Dphip = getDphip3DP3N(Dlambda,lambda(p,:));
       case 'P3jia5B4'
           Dphip = getDphipP3jiapos(Dlambda,lambda(p,:));
       case '3DP2'
            Dphip = getDphip3DP2(Dlambda,lambda(p,:));
       case '3DP3'
            Dphip = getDphip3DP3(Dlambda,lambda(p,:));
       case '3DP3jia8B4'
           Dphip = getDphipP3jia8B4(Dlambda,lambda(p,:));
   end

    
       for i = 1:elemsize 
             for j = 1:elemsizep
                 
                   Dxij =  w(p)*Dphip(:,1,i)*phipp(p,j);
                   Dyij =  w(p)*Dphip(:,2,i)*phipp(p,j);
                   Dzij =  w(p)*Dphip(:,3,i)*phipp(p,j);
                 
                   Dx = Dx + sparse(double(elemdofp(:,j)),double(elemdof(:,i)),Dxij.*volume,Np,Nu);
                   Dy = Dy + sparse(double(elemdofp(:,j)),double(elemdof(:,i)),Dyij.*volume,Np,Nu);
                   Dz = Dz + sparse(double(elemdofp(:,j)),double(elemdof(:,i)),Dzij.*volume,Np,Nu);
                   
             end
      end
end
 B = -[Dx Dy Dz];
end


function [f1,f2,f3]=righthand(Nu,option,pde,NT,node,elem,elemdof,volume,elemsize)
%% 计算右端项
f1 = zeros(Nu,1);
f2 = zeros(Nu,1);
f3 = zeros(Nu,1);
if ~isfield(option,'fquadorder')
    switch option.usolver
        case 'P2jia'
            fquadorder = 5;
        case '3DP2N'
            fquadorder = 5;
        case '3DP3N'
            fquadorder = 6;   % default order
        case 'P3jia5B4'
            fquadorder = 6;
        case '3DP2'
            fquadorder = 5;
        case '3DP3'
            fquadorder = 6;
        case '3DP3jia8B4'
            fquadorder = 6;
    end
end
if ~isfield(pde,'f') || (isreal(pde.f) && (pde.f==0))
    pde.f = [];
end
if ~isempty(pde.f) 
    % quadrature points in the barycentric coordinate
    [lambda,weight] = quadpts3(fquadorder);
    % basis values at quadrature points
   switch option.usolver
       case 'P2jia'
           phi = getphiP2jia(lambda);
       case '3DP2N'
           phi = getphi3DP2N(lambda);
       case '3DP3N'
            phi = getphi3DP3N(lambda);
       case 'P3jia5B4'
           phi = getphiP3jiapos(lambda);  
       case '3DP2'
            phi = getphi3DP2(lambda);
       case '3DP3'
            phi = getphi3DP3(lambda);
       case '3DP3jia8B4'
            phi = getphiP3jia8B4(lambda);  
   end
     
     nQuad = size(lambda,1);
    ft1 = zeros(NT,elemsize);
    ft2 = zeros(NT,elemsize);
    ft3 = zeros(NT,elemsize);
    for p = 1:nQuad
        % quadrature points in the x-y-z coordinate
        pxyz =    lambda(p,1)*node(elem(:,1),:) ...
                    + lambda(p,2)*node(elem(:,2),:) ...
                    + lambda(p,3)*node(elem(:,3),:)...
                    + lambda(p,4)*node(elem(:,4),:);
        % function values at quadrature points
        fp = pde.f(pxyz);
        % evaluate fp outside.
        for j = 1:elemsize
            ft1(:,j) = ft1(:,j) + fp(:,1).*phi(p,j)*weight(p);
            ft2(:,j) = ft2(:,j) + fp(:,2).*phi(p,j)*weight(p);
            ft3(:,j) = ft3(:,j) + fp(:,3).*phi(p,j)*weight(p);
        end
    end
    ft1 = ft1.*repmat(volume,1,elemsize);
    ft2 = ft2.*repmat(volume,1,elemsize);
    ft3 = ft3.*repmat(volume,1,elemsize);
    f1 = accumarray(elemdof(:),ft1(:),[Nu 1]);
    f2 = accumarray(elemdof(:),ft2(:),[Nu 1]);
    f3 = accumarray(elemdof(:),ft3(:),[Nu 1]);
end
end



function [AD,BD,f,u,po,g,ufreeDof,pDof] = getbdStokes(node,eleminforming,bdFlag,f1,f2,f3,...
                                                                                      Np,Nu,elemdof,A,B,pde,elem,option)
%% Boundary condition of Stokes equation

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
        switch option.usolver
            case 'P2jia'
                isFixedDof(elemdof(bdFlag(:,1) == 1,[1,2,3])) = true;
                isFixedDof(elemdof(bdFlag(:,2) == 1,[4,5,6])) = true;        
                isFixedDof(elemdof(bdFlag(:,3) == 1,[7,8,9])) = true;
                isFixedDof(elemdof(bdFlag(:,4) == 1,[10,11,12])) = true;
            case '3DP2N'
                isFixedDof(elemdof(bdFlag(:,1) == 1,[2,3,4,8,9,10])) = true;
                isFixedDof(elemdof(bdFlag(:,2) == 1,[1,3,4,6,7,10])) = true;
                isFixedDof(elemdof(bdFlag(:,3) == 1,[1,2,4,5,7,9])) = true;
                isFixedDof(elemdof(bdFlag(:,4) == 1,[1,2,3,5,6,8])) = true;
                
            case '3DP3N'
                isFixedDof(elemdof(bdFlag(:,1) == 1,[2,3,4,5,15,16,17,18,19,20])) = true;    
                isFixedDof(elemdof(bdFlag(:,2) == 1,[1,3,4,6,11,12,13,14,19,20])) = true;  
                isFixedDof(elemdof(bdFlag(:,3) == 1,[1,2,4,7,9,10,13,14,17,18])) = true;        
                isFixedDof(elemdof(bdFlag(:,4) == 1,[1,2,3,8,9,10,11,12,15,16])) = true;
                
            case 'P3jia5B4'
                isFixedDof(elemdof(bdFlag(:,1) == 1,1:6)) = true;
                isFixedDof(elemdof(bdFlag(:,2) == 1,7:12)) = true;
                isFixedDof(elemdof(bdFlag(:,3) == 1,13:18)) = true;
                isFixedDof(elemdof(bdFlag(:,4) == 1,19:24)) = true;
                
            case '3DP2'
                isFixedDof(elemdof(bdFlag(:,1) == 1,[2,3,4,8,9,10])) = true;
                isFixedDof(elemdof(bdFlag(:,2) == 1,[1,3,4,6,7,10])) = true;
                isFixedDof(elemdof(bdFlag(:,3) == 1,[1,2,4,5,7,9])) = true;
                isFixedDof(elemdof(bdFlag(:,4) == 1,[1,2,3,5,6,8])) = true;
                
            case '3DP3'
                isFixedDof(elemdof(bdFlag(:,1) == 1,[2,3,4,5,15,16,17,18,19,20])) = true;    
                isFixedDof(elemdof(bdFlag(:,2) == 1,[1,3,4,6,11,12,13,14,19,20])) = true;  
                isFixedDof(elemdof(bdFlag(:,3) == 1,[1,2,4,7,9,10,13,14,17,18])) = true;        
                isFixedDof(elemdof(bdFlag(:,4) == 1,[1,2,3,8,9,10,11,12,15,16])) = true;   
                
            case '3DP3jia8B4'
                 isFixedDof(elemdof(bdFlag(:,1) == 1,1:6)) = true;
                isFixedDof(elemdof(bdFlag(:,2) == 1,7:12)) = true;
                isFixedDof(elemdof(bdFlag(:,3) == 1,13:18)) = true;
                isFixedDof(elemdof(bdFlag(:,4) == 1,19:24)) = true;

                
        end
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
        switch option.usolver
            case 'P2jia'
                [lambda,w] = quadpts(5);
            case '3DP2N'
                [lambda,w] = quadpts(5);
            case '3DP3N'
                [lambda,w] = quadpts(6);
            case 'P3jia5B4'
                [lambda,w] = quadpts(6);
            case '3DP2'
                [lambda,w] = quadpts(5);
            case '3DP3'
                [lambda,w] = quadpts(6);
            case '3DP3jia8B4'
                [lambda,w] = quadpts(6);
        end
        nQuad = size(lambda,1);
        % quadratic 

        totalface(:,:,1) = elem(:,[2,3,4]);
        totalface(:,:,2) = elem(:,[3,4,1]);
        totalface(:,:,3) = elem(:,[4,1,2]);
        totalface(:,:,4) = elem(:,[1,2,3]);
        
        
       for j = 1:4
           isNeumann = (bdFlag(:,j)==2);
           Neumannidx = elemdof(isNeumann,:);
          
           if ~isempty(Neumannidx)
                 
                    Neumannidx = Neumannidx(:);
                    Neumann = totalface(:,:,j);
                    Neumann = Neumann(isNeumann,:);
                    areafaces = getareaface(node,Neumann);
   
               for pp = 1:nQuad
                   lam = lambda(pp,:);
                   bdphiface = getbdphiface(lam,j,option);

                   pxyz = lam(1)*node(Neumann(:,1),:)+lam(2)*node(Neumann(:,2),:)...
                                      +lam(3)*node(Neumann(:,3),:);
                                     
                   gp = pde.g_N(pxyz);

                   
                   g1 = w(pp)*areafaces.*gp(:,1).*bdphiface;
                   g2 = w(pp)*areafaces.*gp(:,2).*bdphiface;
                   g3 = w(pp)*areafaces.*gp(:,3).*bdphiface;
                   
                   
                        
                                  
                    S = size(f1,1);
                    f1 = f1 + accumarray(Neumannidx,g1(:),[S,1]);
                    f2 = f2 + accumarray(Neumannidx,g2(:),[S,1]);
                    f3 = f3 + accumarray(Neumannidx,g3(:),[S,1]);
      
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
        u1 = zeros(Nu,1);
        u2 = zeros(Nu,1);
        u3 = zeros(Nu,1);  
        switch option.usolver
            
            case 'P2jia'
                face = eleminforming.face;
                NF = length(face);
                idx = (fixedDof<=NF);
                isDirichlet = fixedDof(idx);
                Dirichletface = face(isDirichlet,:);
                
                 [lambda,w] = quadpts(5); %quador point
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
                
            
            case '3DP2N'
                 edge = eleminforming.edge;
                 N = size(node,1);
                 idx = (fixedDof > N);                         % index of edge nodes
                 uD = pde.g_D(node(fixedDof(~idx),:));                
                 u1(fixedDof(~idx)) = uD(:,1);
                 u2(fixedDof(~idx)) = uD(:,2);
                 u3(fixedDof(~idx)) = uD(:,3);
                 
                 bdEdgeIdx = fixedDof(idx) - N;
                 bdEdgeMid = (node(edge(bdEdgeIdx,1),:) + node(edge(bdEdgeIdx,2),:))/2;
                 uD = pde.g_D(bdEdgeMid);
                 u1(fixedDof(idx)) = uD(:,1);
                 u2(fixedDof(idx)) = uD(:,2);
                 u3(fixedDof(idx)) = uD(:,3);

            
                     
            case '3DP3N'
                ND = max(elem(:));
                NF = max(max(elemdof(:,1:8)));

                idx = (fixedDof<ND+1); 
                uD = pde.g_D(node(fixedDof(idx),:));
                
                u1(fixedDof(idx)) = uD(:,1);
                u2(fixedDof(idx)) = uD(:,2);
                u3(fixedDof(idx)) = uD(:,3);
                
                idx = fixedDof>NF;
                bdEdgeIdx = fixedDof(idx)-NF;
                edges = [eleminforming.edge(:,[1,2]);eleminforming.edge(:,[2,1])];
                
                uD = pde.g_D((2*node(edges(bdEdgeIdx,1),:)+node(edges(bdEdgeIdx,2),:))/3);
                
                u1(fixedDof(idx)) = uD(:,1);
                u2(fixedDof(idx)) = uD(:,2);
                u3(fixedDof(idx)) = uD(:,3);
                
                
            case 'P3jia5B4'
                face = eleminforming.face;
                NF = length(face);
                idx = (fixedDof<=NF);
                isDirichlet = fixedDof(idx);
                Dirichletface = face(isDirichlet,:);
                u1 = zeros(Nu,1);
                u2 = zeros(Nu,1);
                u3 = zeros(Nu,1);    
                
               [lambda,w] = quadpts(6);
               nQuad=length(lambda);
                for p=1:nQuad
                     pxyz=lambda(p,1)*node(Dirichletface(:,1),:)+lambda(p,2)*node(Dirichletface(:,2),:)...
                                        +lambda(p,3)*node(Dirichletface(:,3),:);
                     uD=pde.g_D(pxyz);
                     uD1 = uD(:,1);
                     uD2 = uD(:,2);
                     uD3 = uD(:,3);

                    la1 = lambda(p,1);
                    la2 = lambda(p,2);
                    la3 = lambda(p,3);


                    uD1 = w(p)*[uD1*la1.^2;uD1*la2.^2;uD1*la3.^2;uD1*la1*la2;uD1*la1*la3;uD1*la2*la3];
                    uD2 = w(p)*[uD2*la1.^2;uD2*la2.^2;uD2*la3.^2;uD2*la1*la2;uD2*la1*la3;uD2*la2*la3];
                    uD3 = w(p)*[uD3*la1.^2;uD3*la2.^2;uD3*la3.^2;uD3*la1*la2;uD3*la1*la3;uD3*la2*la3];

                     u1(fixedDof)=u1(fixedDof)+uD1;
                     u2(fixedDof)=u2(fixedDof)+uD2;
                     u3(fixedDof)=u3(fixedDof)+uD3;
                end
                 
                
            case '3DP2'
                 edge = eleminforming.edge;
                 N = size(node,1);
                 idx = (fixedDof > N);                         % index of edge nodes
                 uD = pde.g_D(node(fixedDof(~idx),:));                
                 u1(fixedDof(~idx)) = uD(:,1);
                 u2(fixedDof(~idx)) = uD(:,2);
                 u3(fixedDof(~idx)) = uD(:,3);
                 
                 bdEdgeIdx = fixedDof(idx) - N;
                 bdEdgeMid = (node(edge(bdEdgeIdx,1),:) + node(edge(bdEdgeIdx,2),:))/2;
                 uD = pde.g_D(bdEdgeMid);
                 u1(fixedDof(idx)) = uD(:,1);
                 u2(fixedDof(idx)) = uD(:,2);
                 u3(fixedDof(idx)) = uD(:,3);
                 
                 
                 
            case '3DP3'
                ND = max(elem(:));
                NF = max(max(elemdof(:,1:8)));

                idx = (fixedDof<ND+1); 
                uD = pde.g_D(node(fixedDof(idx),:));
                
                u1(fixedDof(idx)) = uD(:,1);
                u2(fixedDof(idx)) = uD(:,2);
                u3(fixedDof(idx)) = uD(:,3);
                
                idx = fixedDof>NF;
                bdEdgeIdx = fixedDof(idx)-NF;
                edges = [eleminforming.edge(:,[1,2]);eleminforming.edge(:,[2,1])];
                
                uD = pde.g_D((2*node(edges(bdEdgeIdx,1),:)+node(edges(bdEdgeIdx,2),:))/3);
                
                u1(fixedDof(idx)) = uD(:,1);
                u2(fixedDof(idx)) = uD(:,2);
                u3(fixedDof(idx)) = uD(:,3);
                
                
             case '3DP3jia8B4'
                face = eleminforming.face;
                NF = length(face);
                idx = (fixedDof<=NF);
                isDirichlet = fixedDof(idx);
                Dirichletface = face(isDirichlet,:);
                u1 = zeros(Nu,1);
                u2 = zeros(Nu,1);
                u3 = zeros(Nu,1);    
                
               [lambda,w] = quadpts(6);
               nQuad=length(lambda);
                for p=1:nQuad
                     pxyz=lambda(p,1)*node(Dirichletface(:,1),:)+lambda(p,2)*node(Dirichletface(:,2),:)...
                                        +lambda(p,3)*node(Dirichletface(:,3),:);
                     uD=pde.g_D(pxyz);
                     uD1 = uD(:,1);
                     uD2 = uD(:,2);
                     uD3 = uD(:,3);

                    la1 = lambda(p,1);
                    la2 = lambda(p,2);
                    la3 = lambda(p,3);


                    uD1 = w(p)*[uD1*la1.^2;uD1*la2.^2;uD1*la3.^2;uD1*la1*la2;uD1*la1*la3;uD1*la2*la3];
                    uD2 = w(p)*[uD2*la1.^2;uD2*la2.^2;uD2*la3.^2;uD2*la1*la2;uD2*la1*la3;uD2*la2*la3];
                    uD3 = w(p)*[uD3*la1.^2;uD3*la2.^2;uD3*la3.^2;uD3*la1*la2;uD3*la1*la3;uD3*la2*la3];

                     u1(fixedDof)=u1(fixedDof)+uD1;
                     u2(fixedDof)=u2(fixedDof)+uD2;
                     u3(fixedDof)=u3(fixedDof)+uD3;
                end
                
                
                
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





function bdphiface = getbdphiface(lambda,j,option)

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

switch option.usolver
    case 'P2jia'
        bdphiface = getphiP2jia(lam);
        
    case '3DP2N'
        bdphiface = getphi3DP2N(lam);
        
    case '3DP3N'
        bdphiface = getphi3DP3N(lam);
        
    case 'P3jia5B4'
        bdphiface = getphiP3jiapos(lam);
        
    case '3DP2'
        bdphiface = getphi3DP2(lam);
        
    case '3DP3'
        bdphiface = getphi3DP3(lam);
        
    case '3DP3jia8B4'
        bdphiface =  getphiP3jia8B4(lam);
        
end

end


function areaface = getareaface(node,Neumann)
%%计算面的面积


    v12 = node(Neumann(:,2),:)-node(Neumann(:,1),:);
    v13 = node(Neumann(:,3),:)-node(Neumann(:,1),:);
    
    
    areaface = 0.5*sqrt(abs(sum(mycross(v12,v13,2).^2,2)));

end





