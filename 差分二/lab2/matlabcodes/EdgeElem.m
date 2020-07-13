function Edge_elem_idx = EdgeElem(edge)
%% find the elem as point edge for edge
global elem

idx = sum((uint32(elem == edge(1))+uint32(elem == edge(2))),2);
Edge_elem_idx = sort(find(idx==2));