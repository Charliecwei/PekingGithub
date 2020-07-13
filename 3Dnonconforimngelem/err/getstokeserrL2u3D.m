function  err= getstokeserrL2u3D(node,elem,exactu,uh,var)
%%3D stokes u µÄ L2 Îó²î
err = 0;

switch var
    case 'P2jia'   %P2jia3B3 elem 

         elem3dof=dof3P2jia(elem);
         [lambda,w] = quadpts3(4);
         phi = getphiP2jia(lambda);

    case 'P3jia5B4' %P3jia5B4 elem

        elem3dof=dof3P3jia5B4s(elem);
         [lambda,w] = quadpts3(6);
         phi = getphiP3jiapos(lambda);
     
    case '3DP3N'  %P3N elem
         elem3dof=dof3DP3N(elem);
         [lambda,w] = quadpts3(6);
         phi = getphi3DP3N(lambda);
         
    case '3DP2N' %P2N elem
        
        elem3dof=dof3DP2N(elem);
         [lambda,w] = quadpts3(4);
         phi = getphi3DP2N(lambda);
         
    case '3DP2'
        elem3dof=dof3P2(elem);
        [lambda,w] = quadpts3(4);
        phi = getphi3DP2(lambda);
        
    case '3DP3'
        elem3dof=dof3P3(elem);
        [lambda,w] = quadpts3(6);
        phi = getphi3DP3(lambda);
        
     case '3DP3jia8B4' %P3jia5B4 elem

        elem3dof=dof3P3jia8B4s(elem);
        [lambda,w] = quadpts3(6);
        phi = getphiP3jia8B4(lambda);
        
        


end





     nQuad = size(lambda,1);
     Nu=max(elem3dof(:));
     uh1=uh(1:Nu);
     uh2=uh(Nu+1:2*Nu);
     uh3=uh(2*Nu+1:3*Nu);
     
     uh1=uh1(elem3dof);
     uh2=uh2(elem3dof);
     uh3=uh3(elem3dof);
     
     
     for p=1:nQuad
        pxyz=lambda(p,1)*node(elem(:,1),:)+lambda(p,2)*node(elem(:,2),:)...
                                +lambda(p,3)*node(elem(:,3),:)+lambda(p,4)*node(elem(:,4),:);
        u=exactu(pxyz);
        ux=u(:,1);
        uy=u(:,2);
        uz=u(:,3);
        
        
        uhx=uh1*phi(p,:)';
        uhy=uh2*phi(p,:)';
        uhz=uh3*phi(p,:)';
        
        
        err=err+w(p)*((ux-uhx).^2+(uy-uhy).^2+(uz-uhz).^2);
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
     