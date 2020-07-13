function Flip(J)
%% flip inner edge J
global elem iner_edge iner_edge_elem_idx node_iner_edge_idx node_elem_idx

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

node_elem_idx{ad_node(1)} = setdiff(node_elem_idx{ad_node(1)},elem_idx(2));
node_elem_idx{ad_node(2)} = setdiff(node_elem_idx{ad_node(2)},elem_idx(1));


node_elem_idx{ap_node(1)} = [node_elem_idx{ap_node(1)},elem_idx(2)];
node_elem_idx{ap_node(2)} = [node_elem_idx{ap_node(2)},elem_idx(1)];

idx = find(uint32(iner_edge(:,1)==min(ad_node(1),ap_node(2)))+uint32(iner_edge(:,2)==max(ad_node(1),ap_node(2)))==2);
iner_edge_elem_idx(idx,:) = [setdiff(iner_edge_elem_idx(idx,:),elem_idx(2)),elem_idx(1)];

idx = find(uint32(iner_edge(:,1)==min(ad_node(2),ap_node(1)))+uint32(iner_edge(:,2)==max(ad_node(2),ap_node(1)))==2);
iner_edge_elem_idx(idx,:) = [setdiff(iner_edge_elem_idx(idx,:),elem_idx(1)),elem_idx(2)];

