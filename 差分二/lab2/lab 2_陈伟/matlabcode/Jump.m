function jump = Jump(J)
%%find the jump for iner edge J
global node elem iner_edge  iner_edge_elem_idx u0 ;

Ne = length(J);


Con_elem_idx = iner_edge_elem_idx(J,:);

Edge = iner_edge(J,:);

%% T plus
Elem = elem(Con_elem_idx(:,1),:);
Dphip = gradbasis(node,Elem);

Du_p =  u0(Elem(:,1)).*Dphip(:,:,1) + u0(Elem(:,2)).*Dphip(:,:,2)...
              + u0(Elem(:,3)).*Dphip(:,:,3);
 

          
n_p = zeros(Ne,2);
for i = 1:Ne
     idx = find(~((Elem(i,:) == Edge(i,1)) + (Elem(i,:) == Edge(i,2))));
     n_p(i,:) = Dphip(i,:,idx)/(norm(Dphip(i,:,idx),2));
end




%% T minus
Elem = elem(Con_elem_idx(:,2),:);
Dphip = gradbasis(node,Elem);

Du_n =  u0(Elem(:,1)).*Dphip(:,:,1) + u0(Elem(:,2)).*Dphip(:,:,2)...
              + u0(Elem(:,3)).*Dphip(:,:,3);
 

          
n_n = zeros(Ne,2);
for i = 1:Ne
     idx = find(~((Elem(i,:) == Edge(i,1)) + (Elem(i,:) == Edge(i,2))));
     n_n(i,:) = Dphip(i,:,idx)/(norm(Dphip(i,:,idx),2));
end


%% jump
jump = sum(Du_p.*n_p,2) + sum(Du_n.*n_n,2);









