function Inducled_meshs(I)
%% Induced the bound node I  mesh to convex

global  node elem node_iner_edge_idx iner_edge N0;

  thresold =Thresold(I);
  idx = find(thresold<0);
   while ~isempty(idx)
       Edge_idx = node_iner_edge_idx{I}(idx);
       Edge = iner_edge(Edge_idx,:);
       idx = find(uint32(Edge(:,1)>N0)+uint32(Edge(:,2)>N0)==2);                  
%       showmesh(node,elem);
       if ~isempty(idx)
           Flip(Edge_idx(idx(1)));
       else
            Flip(Edge_idx(1));
       end
%        showmesh(node,elem);
        thresold =Thresold(I);
         idx = find(thresold<0);
   end
        
        


