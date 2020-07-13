function delta = Find_deltas(I,thresoldH,b)
global node elem node_elem_idx u0  ;

%% find the delta that u0(I) can ascension and it's convhull area
nelem =  node_elem_idx{I};
Elem = elem(nelem,:);

Dphip = gradbasis(node,Elem);



idx = (Elem == I);


 Dphips(:,:,1) =  u0(Elem(:,1)).*Dphip(:,:,1);
 Dphips(:,:,2) =  u0(Elem(:,2)).*Dphip(:,:,2);
 Dphips(:,:,3) =  u0(Elem(:,3)).*Dphip(:,:,3);
 
 Du_hat= Dphips(:,:,1)+Dphips(:,:,2)+Dphips(:,:,3);
 Dlam = zeros(size(Elem,1),2);
 Dlam(idx(:,1),:) = Dphip(idx(:,1),:,1);
 Dlam(idx(:,2),:) = Dphip(idx(:,2),:,2);
 Dlam(idx(:,3),:) = Dphip(idx(:,3),:,3);
 
 
  
 
 
 
% 
% Elem = elem(nelem,:)';
% edge_ap = sort(reshape(Elem(Elem~=I),2,[]),1)';
% 
% ev = (node(edge_ap(:,1),:)+node(edge_ap(:,2),:))-2*node(I,:);
% k = convhull(ev(:,1),ev(:,2));
%  

 

%       Du_hat = Du_hat(k,:);
%       Dlam = Dlam(k,:);
      Du_hat = [Du_hat;Du_hat(1,:)];
      Dlam = [Dlam;Dlam(1,:)];
      
      Dlampv = [-Dlam(1:end-1,2),Dlam(1:end-1,1)];
      Du_hatpv = [-Du_hat(1:end-1,2),Du_hat(1:end-1,1)];
      
      Dlamm = Dlam(2:end,:);
      Du_hatm = Du_hat(2:end,:);
      
      
      A = 0.5*(sum(dot(Dlampv,Dlamm,2)));
      B = 0.5*(sum(dot(Dlampv,Du_hatm,2))+sum(dot(Du_hatpv,Dlamm,2)));
      C = 0.5*sum(dot(Du_hatpv,Du_hatm,2));

      Delta = roots([A,B,C-b]);
      
      delta = [];
      if isreal(Delta)
          for i = 1:length(Delta)
              if Delta(i)<=thresoldH
                  delta = Delta(i);
              end
          end
          if isempty(delta)
              delta = thresoldH;
          end
      else
          delta = thresoldH;
      end
      

      
 end
 

  
       
     
 
