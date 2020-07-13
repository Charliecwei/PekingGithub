function [soln,eqn] = Poisson3P3(node,elem,bdFlag,pde)


[elem3dof3,face,edge] = dof3P3(elem);

NT = size(elem,1);  
Nu =double( max(elem3dof3(:)));


[Dlambda,volume] = gradbasis3(node,elem);

ii = zeros(210*NT,1); jj = zeros(210*NT,1); 
index = 0;
for i=1:20
    for j=i:20
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
%%%计算矩阵sA
nQuad=size(lambda,1);
sA=zeros(210*NT,nQuad);
Da1=Dlambda(:,:,1);
Da2=Dlambda(:,:,2);
Da3=Dlambda(:,:,3);
Da4=Dlambda(:,:,4);

for pp=1:nQuad
    la1 = lambda(pp,1);
    la2 = lambda(pp,2);
    la3 = lambda(pp,3);
    la4 = lambda(pp,4);
    Dphip = getDphip(Da1,Da2,Da3,Da4,la1,la2,la3,la4);

     index = 0;
    for i = 1:20
        for j = i:20
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
%% 计算右端项
f = zeros(Nu,1);


if ~isempty(pde.f) 
    % quadrature points in the barycentric coordinate
    [lambda,weight] = quadpts3(5);
    % basis values at quadrature points
   
     phi = getphi(lambda);
     
     nQuad = size(lambda,1);
    ft = zeros(NT,20);

    for pp = 1:nQuad
        % quadrature points in the x-y-z coordinate
        pxyz =    lambda(pp,1)*node(elem(:,1),:) ...
                    + lambda(pp,2)*node(elem(:,2),:) ...
                    + lambda(pp,3)*node(elem(:,3),:)...
                    + lambda(pp,4)*node(elem(:,4),:);
        % function values at quadrature points
        fp = pde.f(pxyz);
        % evaluate fp outside.
        for j = 1:20
            ft(:,j) = ft(:,j) + fp(:,1)*phi(pp,j)*weight(pp);
        end
    end
    ft = ft.*repmat(volume,1,20);
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



    % Dirichlet boundary conditions
    if ~isempty(fixedDof) && ~isempty(pde.g_D) && ~(isnumeric(pde.g_D) && (pde.g_D == 0))
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uD = pde.g_D(node(fixedDof,:));  % bd value at vertex dofs.................................
        %找边界面上的基函数,并对其赋值,分别是顶点，面,边
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
