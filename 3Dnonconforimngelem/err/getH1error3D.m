function err=getH1error3D(node,elem,exactDu,uh,var)
% µÃµ½ H1°ë·¶Îó²î

err = 0;
if size(node,2)==3
    [Dlambda,volume] = gradbasis3(node,elem);
    
    switch var
 
        case 'P2jia'


            elemdof = dof3P2jia(elem);
            [lambda,w] = quadpts3(6);
            nQuad = length(lambda);
            for p=1:nQuad
                la1 = lambda(p,1);
                la2 = lambda(p,2);
                la3 = lambda(p,3);
                la4 = lambda(p,4);

                Da1 = Dlambda(:,:,1);
                Da2 = Dlambda(:,:,2);
                Da3 = Dlambda(:,:,3);
                Da4 = Dlambda(:,:,4);

                Dphip = getDphipP2jia(Da1,Da2,Da3,Da4,la1,la2,la3,la4);

                Duh=0;

                for j=1:13
                    Duh = Duh+repmat(uh(elemdof(:,j)),1,3).*Dphip(:,:,j);
                end

                pxyz=la1*node(elem(:,1),:)+la2*node(elem(:,2),:)...
                                      +la3*node(elem(:,3),:)+la4*node(elem(:,4),:);


                Du = exactDu(pxyz);  

                err=err+w(p)*sum((Duh-Du).^2,2);
            end


        case 'P3jia5B4s'


            elem3dof = dof3P3jia5B4s(elem);
            [lambda,weight] = quadpts3(6);
            nQuad = size(lambda,1);
            for p = 1:nQuad
                % quadrature points in the x-y-z coordinate
                pxyz = lambda(p,1)*node(elem(:,1),:) ...
                     + lambda(p,2)*node(elem(:,2),:) ...
                     + lambda(p,3)*node(elem(:,3),:) ...
                     + lambda(p,4)*node(elem(:,4),:);



                Dphip = getDphipP3jiapos(Dlambda,lambda(p,:));



                   Duh=0;
                   for i=1:25
                         Duh =Duh+ repmat(uh(elem3dof(:,i)),1,3).*Dphip(:,:,i);
                   end



                   Du = exactDu(pxyz);
                   err = err + weight(p)*sum(((Du-Duh).^2),2);
            end
            
            
            
            case 'P3jia8B4s'


            elem3dof = dof3P3jia8B4s(elem);
            [lambda,weight] = quadpts3(6);
            nQuad = size(lambda,1);
            for p = 1:nQuad
                % quadrature points in the x-y-z coordinate
                pxyz = lambda(p,1)*node(elem(:,1),:) ...
                     + lambda(p,2)*node(elem(:,2),:) ...
                     + lambda(p,3)*node(elem(:,3),:) ...
                     + lambda(p,4)*node(elem(:,4),:);



                Dphip = getDphipP3jia8B4(Dlambda,lambda(p,:));



                   Duh=0;
                   for i=1:28
                         Duh =Duh+ repmat(uh(elem3dof(:,i)),1,3).*Dphip(:,:,i);
                   end



                   Du = exactDu(pxyz);
                   err = err + weight(p)*sum(((Du-Duh).^2),2);
            end



         case 'P2' %P2 elem

            elem3dof = dof3P2(elem);
            [lambda,weight] = quadpts3(4);
            nQuad = size(lambda,1);
            for p = 1:nQuad
                % quadrature points in the x-y-z coordinate
                pxyz = lambda(p,1)*node(elem(:,1),:) ...
                     + lambda(p,2)*node(elem(:,2),:) ...
                     + lambda(p,3)*node(elem(:,3),:) ...
                     + lambda(p,4)*node(elem(:,4),:);



                Dphip = getDphi2(Dlambda,lambda(p,:));



                   Duh=0;
                   for i=1:10
                         Duh =Duh+ repmat(uh(elem3dof(:,i)),1,3).*Dphip(:,:,i);
                   end



                   Du = exactDu(pxyz);
                   err = err + weight(p)*sum(((Du-Duh).^2),2);
            end

        case 'P1' %P1 elem

            elem3dof = elem;
            [lambda,weight] = quadpts3(1);
            nQuad = size(lambda,1);
            for p = 1:nQuad
                % quadrature points in the x-y-z coordinate
                pxyz = lambda(p,1)*node(elem(:,1),:) ...
                     + lambda(p,2)*node(elem(:,2),:) ...
                     + lambda(p,3)*node(elem(:,3),:) ...
                     + lambda(p,4)*node(elem(:,4),:);



                Dphip = Dlambda;



                   Duh=0;
                   for i=1:4
                         Duh =Duh+ repmat(uh(elem3dof(:,i)),1,3).*Dphip(:,:,i);
                   end



                   Du = exactDu(pxyz);
                   err = err + weight(p)*sum(((Du-Duh).^2),2);
            end
            
            
        case '3DP2N'
            elem3dof = dof3DP2N(elem);
             [lambda,weight] = quadpts3(4);
            nQuad = size(lambda,1);
             for p = 1:nQuad
                % quadrature points in the x-y-z coordinate
                pxyz = lambda(p,1)*node(elem(:,1),:) ...
                     + lambda(p,2)*node(elem(:,2),:) ...
                     + lambda(p,3)*node(elem(:,3),:) ...
                     + lambda(p,4)*node(elem(:,4),:);



                Dphip = getDphip3DP2N(Dlambda,lambda(p,:));



                   Duh=0;
                   for i=1:11
                         Duh =Duh+ repmat(uh(elem3dof(:,i)),1,3).*Dphip(:,:,i);
                   end



                   Du = exactDu(pxyz);
                   err = err + weight(p)*sum(((Du-Duh).^2),2);
             end
            
             
        case '3DP3N'
             elem3dof = dof3DP3N(elem);
             [lambda,weight] = quadpts3(5);
            nQuad = size(lambda,1);
             for p = 1:nQuad
                % quadrature points in the x-y-z coordinate
                pxyz = lambda(p,1)*node(elem(:,1),:) ...
                     + lambda(p,2)*node(elem(:,2),:) ...
                     + lambda(p,3)*node(elem(:,3),:) ...
                     + lambda(p,4)*node(elem(:,4),:);


                 
                Dphip = getDphip3DP3N(Dlambda,lambda(p,:));



                   Duh=0;
                   for i=1:21
                         Duh =Duh+ repmat(uh(elem3dof(:,i)),1,3).*Dphip(:,:,i);
                   end



                   Du = exactDu(pxyz);
                   err = err + weight(p)*sum(((Du-Duh).^2),2);
             end
            
            
            
    
    end
    
    
    
elseif size(node,2)==2
    [Dlambda,volume] = gradbasis(node,elem);
    switch var
        case '2DP3N' %2D P3 Nonconforming element
              elem2Ddof = dof2DP3N(elem);
              [lambda,weight] = quadpts(6);
               nQuad = size(lambda,1);
        for p = 1:nQuad
            % quadrature points in the x-y-z coordinate
            pxyz = lambda(p,1)*node(elem(:,1),:) ...
                 + lambda(p,2)*node(elem(:,2),:) ...
                 + lambda(p,3)*node(elem(:,3),:) ;


            la1=lambda(p,1);
            la2=lambda(p,2);
            la3=lambda(p,3);

            Dla1=Dlambda(:,:,1);
            Dla2=Dlambda(:,:,2);
            Dla3=Dlambda(:,:,3);

            Dphip = getDphip2DP3N(Dla1,Dla2,Dla3,la1,la2,la3);



               Duh=0;
               for i=1:10
                     Duh =Duh+ repmat(uh(elem2Ddof(:,i)),1,2).*Dphip(:,:,i);
               end



               Du = exactDu(pxyz);
               err = err + weight(p)*sum(((Du-Duh).^2),2);
        end
        
        case '2DP2N'
                  elem2Ddof = dof2DP2N(elem);
                  [lambda,weight] = quadpts(4);
                   nQuad = size(lambda,1);
            for p = 1:nQuad
                % quadrature points in the x-y-z coordinate
                pxyz = lambda(p,1)*node(elem(:,1),:) ...
                     + lambda(p,2)*node(elem(:,2),:) ...
                     + lambda(p,3)*node(elem(:,3),:) ;


                Dphip = getDphip2DP2N(Dlambda,lambda(p,:));



                   Duh=0;
                   for i=1:7
                         Duh =Duh+ repmat(uh(elem2Ddof(:,i)),1,2).*Dphip(:,:,i);
                   end



                   Du = exactDu(pxyz);
                   err = err + weight(p)*sum(((Du-Duh).^2),2);
            end
            
    end
           
          
end
    
err = volume.*err;
err(isnan(err)) = 0; % singular values are excluded
err = sqrt(sum(err));

end




