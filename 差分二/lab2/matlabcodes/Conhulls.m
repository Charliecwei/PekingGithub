function s = Conhulls(I)
%% Compute the Du on elem for I point and compute it's hull area
global node elem node_elem_idx u0

Elem = elem(node_elem_idx{I},:);
Dphip = gradbasis(node,Elem);

 Dphips(:,:,1) =  u0(Elem(:,1)).*Dphip(:,:,1);
 Dphips(:,:,2) =  u0(Elem(:,2)).*Dphip(:,:,2);
 Dphips(:,:,3) =  u0(Elem(:,3)).*Dphip(:,:,3);
    
 Du= Dphips(:,:,1)+Dphips(:,:,2)+Dphips(:,:,3);

 


Du = [Du;Du(1,:)];


Dupv = Du(1:end-1,:);
Dupv = [-Dupv(:,2),Dupv(:,1)];

s = 0.5*sum(dot(Dupv,Du(2:end,:),2));
 
 


