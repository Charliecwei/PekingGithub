function [soln,eqn] = Poisson3DP3N(node,elem,bdFlag,pde)


[elem3dof3,face,edge] = dof3DP3N(elem);

NT = size(elem,1);  
Nu =double( max(elem3dof3(:)));


[Dlambda,volume] = gradbasis3(node,elem);

ii = zeros(231*NT,1); jj = zeros(231*NT,1); 
index = 0;
for i=1:21
    for j=i:21
         ii(index+1:index+NT) = double(elem3dof3(:,i)); 
         jj(index+1:index+NT) = double(elem3dof3(:,j));  
        index = index + NT;
    end
end


[lambda, w] = quadpts3(5);

sA=Dphips(lambda,Dlambda,w,volume,NT);


diagIdx = (ii == jj);   upperIdx = ~diagIdx;
A = sparse(ii(diagIdx),jj(diagIdx),sA(diagIdx),Nu,Nu);
AU = sparse(ii(upperIdx),jj(upperIdx),sA(upperIdx),Nu,Nu);
A = A + AU + AU';
clear Aij ii jj sA

f=righthand(Nu,pde,NT,node,elem,elem3dof3,volume);

 [AD,f,u,freeDof] = getbdStokesP3jiaP2(node,face,edge,bdFlag,f,...
                                                                                     Nu,elem3dof3,A,elem,pde);


u(freeDof) = AD(freeDof,freeDof)\f(freeDof);


 soln = struct('u',u);
eqn = struct('A',AD,'f',f,'face',face,'freeDof',freeDof,'Lap',A,'edge',edge);


end



function   sA=Dphips(lambda,Dlambda,w,volume,NT)
%%%ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?sA
nQuad=size(lambda,1);
sA=zeros(210*NT,nQuad);
% Da1=Dlambda(:,:,1);
% Da2=Dlambda(:,:,2);
% Da3=Dlambda(:,:,3);
% Da4=Dlambda(:,:,4);

for pp=1:nQuad
%     la1 = lambda(pp,1);
%     la2 = lambda(pp,2);
%     la3 = lambda(pp,3);
%     la4 = lambda(pp,4);
    Dphip = getDphip3DP3N(Dlambda,lambda(pp,:));

     index = 0;
    for i = 1:21
        for j = i:21
            Aij = 0;

                Aij = Aij + w(pp)*dot(Dphip(:,:,i),Dphip(:,:,j),2);

            Aij = Aij.*volume;
            sA(index+1:index+NT,pp) = Aij;
            index = index + NT;
        end
    end
end
     sA = sum(sA,2);
end




function f=righthand(Nu,pde,NT,node,elem,elem3dof3,volume)
%% ï¿½ï¿½ï¿½ï¿½ï¿½Ò¶ï¿½ï¿½ï¿½
f = zeros(Nu,1);


if ~isempty(pde.f) 
    % quadrature points in the barycentric coordinate
    [lambda,weight] = quadpts3(5);
    % basis values at quadrature points
   
     phi = getphi3DP3N(lambda);
     
     nQuad = size(lambda,1);
    ft = zeros(NT,21);

    for pp = 1:nQuad
        % quadrature points in the x-y-z coordinate
        pxyz =    lambda(pp,1)*node(elem(:,1),:) ...
                    + lambda(pp,2)*node(elem(:,2),:) ...
                    + lambda(pp,3)*node(elem(:,3),:)...
                    + lambda(pp,4)*node(elem(:,4),:);
        % function values at quadrature points
        fp = pde.f(pxyz);
        % evaluate fp outside.
        for j = 1:21
            ft(:,j) = ft(:,j) + fp(:,1)*phi(pp,j)*weight(pp);
        end
    end
    ft = ft.*repmat(volume,1,21);
    f = accumarray(elem3dof3(:),ft(:),[Nu 1]);

end
end





function [AD,f,u,ufreeDof] = getbdStokesP3jiaP2(node,face,edge,bdFlag,f,...
                                                                                      Nu,elem3dof3,A,elem,pde)
    %% Boundary condition of Stokes equation: P2jia-P1 elements

    %% Initial set up
    u = zeros(Nu,1);    
   

    %% Part 1: Find Dirichlet dof and modify the matrix
    % Find Dirichlet boundary dof: fixedDof and pDof
    isFixedDof = false(Nu,1);     
    if ~isempty(bdFlag)       % case: bdFlag is not empty 
        isFixedDof(elem3dof3(bdFlag(:,1) == 1,[2,3,4,5,15,16,17,18,19,20])) = true;    
        isFixedDof(elem3dof3(bdFlag(:,2) == 1,[1,3,4,6,11,12,13,14,19,20])) = true;  
        isFixedDof(elem3dof3(bdFlag(:,3) == 1,[1,2,4,7,9,10,13,14,17,18])) = true;        
        isFixedDof(elem3dof3(bdFlag(:,4) == 1,[1,2,3,8,9,10,11,12,15,16])) = true;   
       
        fixedDof = find(isFixedDof);
        ufreeDof = find(~isFixedDof);
    end
 
  
    % Modify the matrix
    % Build Dirichlet boundary condition into the matrix AD by enforcing
    % AD(fixedDof,fixedDof)=I, AD(fixedDof,ufreeDof)=0, AD(ufreeDof,fixedDof)=0.
    % BD(:,fixedDof) = 0 and thus BD'(fixedDof,:) = 0.
    bdidx = zeros(Nu,1); 
    bdidx(fixedDof) = 1;
    Tbd = spdiags(bdidx,0,Nu,Nu);
    T = spdiags(1-bdidx,0,Nu,Nu);
    AD = T*A*T + Tbd;

    Ndof = Nu;


    %% Part 2: Find boundary faces and modify the load b
    % Find boundary faces bdFace for Neumann boundary condition
    bdFace2dof = [];
    if isempty(bdFlag) && (~isempty(pde.g_N) || ~isempty(pde.g_R))
        % no bdFlag, only pde.g_N or pde.g_R is given in the input
        bdFlag = setboundary3(node,elem,'Neumann');        
    end
    if ~isempty(bdFlag)
        % Find boundary edges and nodes
        bdFace2dof = [elem3dof3(bdFlag(:,1) == 2,[2,3,4,5,15,16,17,18,19,20,21]); ...  
                      elem3dof3(bdFlag(:,2) == 2,[3,4,1,6,19,20,12,11,14,13,21]); ...
                      elem3dof3(bdFlag(:,3) == 2,[4,1,2,7,14,13,18,17,9,10,21]); ...
                      elem3dof3(bdFlag(:,4) == 2,[1,2,3,8,9,10,11,12,15,16,21])];
    end
    % Neumann boundary condition
    if ~isempty(bdFace2dof) && ~isempty(pde.g_N) && ~(isnumeric(pde.g_N) && (pde.g_N == 0))
        v12 = node(bdFace2dof(:,2),:)-node(bdFace2dof(:,1),:);
        v13 = node(bdFace2dof(:,3),:)-node(bdFace2dof(:,1),:);
        area = 0.5*sqrt(abs(sum(mycross(v12,v13,2).^2,2)));
        [lambdagN,weightgN] = quadpts(4);
 
        lagN1 = lambdagN(:,1);
        lagN2 = lambdagN(:,2);
        lagN3 = lambdagN(:,3);
        phigN(:,1) = 0.5*lagN1.*(3*lagN1-1).*(3*lagN1-2);
        phigN(:,2) = 0.5*lagN2.*(3*lagN2-1).*(3*lagN2-2);
        phigN(:,3) = 0.5*lagN3.*(3*lagN3-1).*(3*lagN3-2);
        phigN(:,4) = 27*lagN1.*lagN2.*lagN3;

        phigN(:,5) = 4.5*lagN1.*lagN2.*(3*lagN1-1);
        phigN(:,6) = 4.5*lagN2.*lagN1.*(3*lagN2-1);

        phigN(:,7) = 4.5*lagN1.*lagN3.*(3*lagN1-1);
        phigN(:,8) = 4.5*lagN3.*lagN1.*(3*lagN3-1);

        phigN(:,9) = 4.5*lagN2.*lagN3.*(3*lagN2-1);
        phigN(:,10) = 4.5*lagN3.*lagN2.*(3*lagN3-1);

        phigN(:,11) = 9* (lagN1.^2+lagN2.^2+lagN3.^2)...
        -11*(lagN1.^3+lagN2.^3+lagN3.^3)...
        -72*lagN1.*lagN2.*lagN3;



        nQuadgN = size(lambdagN,1);
        gf = zeros(size(bdFace2dof,1),11);
        for pp = 1:nQuadgN
            % quadrature points in the x-y coordinate
            ppxyz = lambdagN(pp,1)*node(bdFace2dof(:,1),:) ...
                  + lambdagN(pp,2)*node(bdFace2dof(:,2),:) ...
                  + lambdagN(pp,3)*node(bdFace2dof(:,3),:);
            gNp = pde.g_N(ppxyz);
            for iN = 1:11
                gf(:,iN) = gf(:,iN) + weightgN(pp)*phigN(pp,iN)*gNp;
            end
        end
        gf = gf.*repmat(area,1,11);
        f = f + accumarray(bdFace2dof(:),gf(:),[Ndof,1]); 
    end









    % Dirichlet boundary conditions
    if ~isempty(fixedDof) && ~isempty(pde.g_D) && ~(isnumeric(pde.g_D) && (pde.g_D == 0))
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uD = pde.g_D(node(fixedDof,:));  % bd value at vertex dofs.................................
        ND = max(elem(:));
        NF = max(max(elem3dof3(:,1:8)));
        
        idx = (fixedDof<ND+1); 
        uD = pde.g_D(node(fixedDof(idx),:));
        
       u(fixedDof(idx)) = uD;

       
       
       idx = and(fixedDof>ND,fixedDof<NF+1);
       bdFaceIdx = fixedDof(idx)-ND;
       uD = pde.g_D((node(face(bdFaceIdx,1),:)+node(face(bdFaceIdx,2),:)+node(face(bdFaceIdx,3),:))/3);
       
       u(fixedDof(idx)) = uD;

        
       
       idx = fixedDof>NF;
       bdEdgeIdx = fixedDof(idx)-NF;
      edges = [edge(:,[1,2]);edge(:,[2,1])];
       
       
       
       uD = pde.g_D((2*node(edges(bdEdgeIdx,1),:)+node(edges(bdEdgeIdx,2),:))/3);
             
             
       u(fixedDof(idx)) = uD;

        
        
        
        
        

        f = f - A*u;  % bring affect of nonhomgenous Dirichlet bd condition to
        
        f(fixedDof) = u(fixedDof);
 
    end
    % The case non-empty Dirichlet but g_D=[] corresponds to the zero Dirichlet
    % boundary condition and no modification is needed.
    end 



