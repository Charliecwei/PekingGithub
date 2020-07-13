function err=getL2error3D(node,elem,exactu,uh,var)
%得到L2误差估计

err = 0;
if size(node,2)==3 % 3D 
    switch var


        case 'P3jia5B4s'  % P3jia5B4s elem
               elemdof = dof3P3jia5B4s(elem);
              %8阶数值积分
              quadOrder = 8;
              [lambda,weight] = quadpts3(quadOrder);
              phi = getphiP3jiapos(lambda);
              
        case 'P3jia8B4s'  % P3jia5B4s elem
               elemdof = dof3P3jia8B4s(elem);
              %8阶数值积分
              quadOrder = 8;
              [lambda,weight] = quadpts3(quadOrder);
              phi = getphiP3jia8B4(lambda);


        case 'P2jia' %P2jia3B3 elem
             elemdof = dof3P2jia(elem);
             %6阶数值积分
              quadOrder = 6;
              [lambda,weight] = quadpts3(quadOrder);
              phi = getphiP2jia(lambda);

        case 'P2'  %P2 elem
             elemdof = dof3P2(elem); 
             quadOrder = 4;
              [lambda,weight] = quadpts3(quadOrder);
              phi = getphi3DP2(lambda);


        case 'P1' %P1 elem
            elemdof = elem;
            quadOrder = 2;
            [lambda,weight] = quadpts3(quadOrder);
             phi = lambda;

        case '3DP2N' %P2 nonconforming elem
            elemdof = dof3DP2N(elem);
             quadOrder = 4;
              [lambda,weight] = quadpts3(quadOrder);
              phi = getphi3DP2N(lambda);
              
        case '3DP3N' %P3 nonconforming elem
             elemdof = dof3DP3N(elem);
             quadOrder = 5;
             [lambda,weight] = quadpts3(quadOrder);
              phi = getphi3DP3N(lambda);
            


    end




     nQuad = size(lambda,1);
    for p = 1:nQuad
        % evaluate uh at quadrature point
         uhp = sum(uh(elemdof(:,:)).*phi(p,:) ,2);
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
    
elseif size(node,2)==2
    
    switch var
        
          case '2DP3N'  % 2D P3 nonconforming
           elemdof = dof2DP3N(elem);
          quadOrder = 6;
          [lambda,weight] = quadpts(quadOrder);
          phi = getphi2DP3N(lambda);
          
        case '2DP2N' %2D  P2 nonconforming
           elemdof = dof2DP2N(elem);
          quadOrder = 4;
          [lambda,weight] = quadpts(quadOrder);
          phi = getphi2DP2N(lambda);
            
          
    end
    
    nQuad = size(lambda,1);
    for p = 1:nQuad
        % evaluate uh at quadrature point
         uhp = sum(uh(elemdof(:,:)).*phi(p,:) ,2);
         % quadrature points in the x-y coordinate
        pxy = lambda(p,1)*node(elem(:,1),:) ...
            + lambda(p,2)*node(elem(:,2),:) ...
            + lambda(p,3)*node(elem(:,3),:) ;
        err = err + weight(p)*((exactu(pxy) - uhp).^2);
    end
    
     %% Modification
    % volume of trangle
   [~,volume] = gradbasis(node,elem);
    err(isnan(err)) = 0; % singular point is excluded
    err = sqrt(sum(volume.*err));
    
          

        

                 



end             



