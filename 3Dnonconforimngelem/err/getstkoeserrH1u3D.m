function err=getstkoeserrH1u3D(node,elem,exactDu,uh,var)
%%3D stokes u µÄ H1 °ë·¶Îó²î
err = 0;
switch var

case 'P2jia' %P2jia3B3 elem

    elem3dof = dof3P2jia(elem);
    [lambda,w] = quadpts3(4);


case 'P3jia5B4'

   elem3dof = dof3P3jia5B4s(elem);
    [lambda,w] = quadpts3(5);



case '3DP3N'
    elem3dof = dof3DP3N(elem);
    [lambda,w] = quadpts3(5);
    
    
case '3DP2N'
    elem3dof = dof3DP2N(elem);
    [lambda,w] = quadpts3(4);
    
case '3DP2'
    elem3dof = dof3P2(elem);
    [lambda,w] = quadpts3(4);
    
case '3DP3' 
    elem3dof = dof3P3(elem);
    [lambda,w] = quadpts3(5);
        
   
case '3DP3jia8B4'

   elem3dof = dof3P3jia8B4s(elem);
    [lambda,w] = quadpts3(5);


       

    
    
end

    Nu = max(elem3dof(:));
    uh1 = uh(1:Nu);
    uh2 = uh(Nu+1:2*Nu);
    uh3 = uh(2*Nu+1:3*Nu);

    [Dlambda,volume] = gradbasis3(node,elem);

    nQuad = length(lambda);



    for p=1:nQuad
        
        switch var
            case 'P2jia'

                Dphip = getDphipP2jia(Dlambda,lambda(p,:)); 
                
            case 'P3jia5B4'

                Dphip = getDphipP3jiapos(Dlambda,lambda(p,:));
                
            case '3DP3N'
                Dphip = getDphip3DP3N(Dlambda,lambda(p,:));
                
            case '3DP2N'
                Dphip = getDphip3DP2N(Dlambda,lambda(p,:));
                
            case '3DP2'
                 Dphip = getDphip3DP2(Dlambda,lambda(p,:));
                 
            case '3DP3'
                 Dphip = getDphip3DP3(Dlambda,lambda(p,:));
                 
            case '3DP3jia8B4'
                Dphip = getDphipP3jia8B4(Dlambda,lambda(p,:));
                 
                   
        end

        Duh1=0;
        Duh2=0;
        Duh3=0;

        for j=1:size(elem3dof,2)
            Duh1 = Duh1+repmat(uh1(elem3dof(:,j)),1,3).*Dphip(:,:,j);
            Duh2 = Duh2+repmat(uh2(elem3dof(:,j)),1,3).*Dphip(:,:,j);
            Duh3 = Duh3+repmat(uh3(elem3dof(:,j)),1,3).*Dphip(:,:,j);
        end

        pxyz=lambda(p,1)*node(elem(:,1),:)+lambda(p,2)*node(elem(:,2),:)...
                              +lambda(p,3)*node(elem(:,3),:)+lambda(p,4)*node(elem(:,4),:);


        Du1 = exactDu.u1(pxyz);  
        Du2 = exactDu.u2(pxyz);
        Du3 = exactDu.u3(pxyz);

        err=err+w(p)*sum((Duh1-Du1).^2+(Duh2-Du2).^2+(Duh3-Du3).^2,2);
    end
    


err = volume.*err;
err(isnan(err)) = 0; % singular values are excluded
err = sqrt(sum(err));



end



