function Flips(J)
%% flip inner edge J
global elem node iner_edge iner_edge_elem_idx node_iner_edge_idx node_elem_idx

ad_node = iner_edge(J,:);
elem_idx = iner_edge_elem_idx(J,:);
elem_edge = elem(elem_idx,:);

ap_node(1) = setdiff(elem_edge(1,:),ad_node);
ap_node(2) = setdiff(elem_edge(2,:),ad_node);

iner_edge(J,:) = sort(ap_node);

elem(elem_idx,:) = [ad_node(1),ap_node;ad_node(2),ap_node];

node_iner_edge_idx{ad_node(1)} = setdiff(node_iner_edge_idx{ad_node(1)},J);
node_iner_edge_idx{ad_node(2)} = setdiff(node_iner_edge_idx{ad_node(2)},J);

node_iner_edge_idx{ap_node(1)} = [node_iner_edge_idx{ap_node(1)},J];
node_iner_edge_idx{ap_node(2)} = [node_iner_edge_idx{ap_node(2)},J];

node_elem_idx{ad_node(1)} = setdiff(node_elem_idx{ad_node(1)},elem_idx(2),'stable');
node_elem_idx{ad_node(2)} = setdiff(node_elem_idx{ad_node(2)},elem_idx(1),'stable');

%%
ev = [node(ad_node(1),:)+node(ap_node(2),:)-2*node(ap_node(1),:);...
          node(ad_node(2),:)+node(ap_node(2),:)-2*node(ap_node(1),:)];
ev(2,:) = [-ev(2,2),ev(2,1)];
s = dot(ev(1,:),ev(2,:),2);
if s>0
    npeidx = node_elem_idx{ap_node(1)};
     idx = (npeidx==elem_idx(1));
     idx_g = find(idx);
     if idx_g==1
         node_elem_idx{ap_node(1)} = [npeidx,elem_idx(2)];
     else
          node_elem_idx{ap_node(1)} = [npeidx(1:idx_g-1),elem_idx(2),npeidx(idx_g:end)];
     end
     
      npeidx = node_elem_idx{ap_node(2)};
      idx = (npeidx==elem_idx(2));
      idx_g = find(idx);
      if idx_g == 1
          node_elem_idx{ap_node(2)} = [npeidx,elem_idx(1)];
      else
          node_elem_idx{ap_node(2)} = [npeidx(1:idx_g-1),elem_idx(1),npeidx(idx_g:end)];
      end
      
else
     npeidx = node_elem_idx{ap_node(1)};
     idx = (npeidx==elem_idx(1));
     idx_g = find(idx);
     if idx_g == length(idx)
         node_elem_idx{ap_node(1)} = [elem_idx(2),npeidx];
     else
         node_elem_idx{ap_node(1)} = [npeidx(1:idx_g),elem_idx(2),npeidx(idx_g+1:end)];
     end
      npeidx = node_elem_idx{ap_node(2)};
      idx = (npeidx==elem_idx(2));
      idx_g = find(idx);
      if idx_g == length(idx)
          node_elem_idx{ap_node(2)} = [elem_idx(1),npeidx];
      else
          node_elem_idx{ap_node(2)} = [npeidx(1:idx_g),elem_idx(1),npeidx(idx_g+1:end)];
      end
      
end
    
          

     
%      
% node_elem_idx{ap_node(1)} = [node_elem_idx{ap_node(1)},elem_idx(2)];
% node_elem_idx{ap_node(2)} = [node_elem_idx{ap_node(2)},elem_idx(1)];

idx = find(uint32(iner_edge(:,1)==min(ad_node(1),ap_node(2)))+uint32(iner_edge(:,2)==max(ad_node(1),ap_node(2)))==2);
iner_edge_elem_idx(idx,:) = [setdiff(iner_edge_elem_idx(idx,:),elem_idx(2)),elem_idx(1)];

idx = find(uint32(iner_edge(:,1)==min(ad_node(2),ap_node(1)))+uint32(iner_edge(:,2)==max(ad_node(2),ap_node(1)))==2);
iner_edge_elem_idx(idx,:) = [setdiff(iner_edge_elem_idx(idx,:),elem_idx(1)),elem_idx(2)];

