function Inducled_mesh(I)
%% Induced the bound node I  mesh to convex

global  elem iner_edge node_elem_idx iner_edge_elem_idx u0 ubd N0;




Elem = elem(node_elem_idx{I},:)';
edge_appose = sort(reshape(Elem(Elem~=I),2,[]),1)';  %% the egde as i for vecter
Elem_idx = [];
for i = 1:size(edge_appose,1)
    Js = uint16(iner_edge(:,1)==edge_appose(i,1))+uint16(iner_edge(:,2)==edge_appose(i,2))==2;
    if ~isempty(iner_edge_elem_idx(Js,:))
         Elem_idx = [Elem_idx,iner_edge_elem_idx(Js,:)];
    end
end
Elem_idx = unique(Elem_idx);

Elem = elem(Elem_idx,:);

totalEdge = uint16([Elem(:,[2,3]); Elem(:,[3,1]); Elem(:,[1,2])]);
totalEdges = sort(totalEdge,2);
[Edge, ~, ~] = unique(totalEdges,'rows','legacy');

J = [];
for i = 1:size(Edge,1)
    Js = find(uint16(iner_edge(:,1)==Edge(i,1))+uint16(iner_edge(:,2)==Edge(i,2)));
    J = [J;Js];
end
J = unique(J);




jump = Jump(J);
idxj = find(jump<0);
 


while ~isempty(idxj)
    
  Js = J(idxj);
   K = isconve(Js); 
    
   if sum(K)==0
       Js = Js(1);
       Edge = iner_edge(Js,:);
       Elems = elem(iner_edge_elem_idx(Js,:),:);

       appose_node = setdiff(Elems,Edge);
       idx = appose_node>I;
       Ij = max(appose_node(idx));

       u0(Ij) = ubd(Ij-N0);
       
   else
       
   
    idx = find(K>0,1);
    Js = Js(idx);
  
    
    Flips(Js);
    
    
  
    
   end
    %%

    
   Elem = elem(node_elem_idx{I},:)';
    edge_appose = sort(reshape(Elem(Elem~=I),2,[]),1)';  %% the egde as i for vecter
    Elem_idx = [];
    for i = 1:size(edge_appose,1)
        Js = uint16(iner_edge(:,1)==edge_appose(i,1))+uint16(iner_edge(:,2)==edge_appose(i,2))==2;
          if ~isempty(iner_edge_elem_idx(Js,:))
             Elem_idx = [Elem_idx,iner_edge_elem_idx(Js,:)];
         end
  end
    Elem_idx = unique(Elem_idx);

    Elem = elem(Elem_idx,:);

    totalEdge = uint16([Elem(:,[2,3]); Elem(:,[3,1]); Elem(:,[1,2])]);
    totalEdges = sort(totalEdge,2);
    [Edge, ~, ~] = unique(totalEdges,'rows','legacy');

    J = [];
    for i = 1:size(Edge,1)
        Js = find(uint16(iner_edge(:,1)==Edge(i,1))+uint16(iner_edge(:,2)==Edge(i,2)));
        J = [J;Js];
    end
    J = unique(J);
    
    %%
    jump = Jump(J);
    idxj = find(jump<0);
 
end


 
end

