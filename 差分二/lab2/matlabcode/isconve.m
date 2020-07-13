function s = isconve(J)
global node elem iner_edge iner_edge_elem_idx
%% judge the tow elem have a com edge is convex 
%%s = 1 convex    s=0, noconvex
%%
s = 0*J;
Ne = length(J);
if ~isempty(J)
    Elem1 = elem(iner_edge_elem_idx(J,1),:);

    %ve1 = node(Elem1(:,3),:)-node(Elem1(:,2),:);
    ve2 = node(Elem1(:,1),:)-node(Elem1(:,3),:);
    ve3 = node(Elem1(:,2),:)-node(Elem1(:,1),:);
    area1 = 0.5*(-ve3(:,1).*ve2(:,2) + ve3(:,2).*ve2(:,1));

    Elem2 = elem(iner_edge_elem_idx(J,2),:);

    %ve1 = node(Elem2(:,3),:)-node(Elem2(:,2),:);
    ve2 = node(Elem2(:,1),:)-node(Elem2(:,3),:);
    ve3 = node(Elem2(:,2),:)-node(Elem2(:,1),:);
    area2 = 0.5*(-ve3(:,1).*ve2(:,2) + ve3(:,2).*ve2(:,1));
    
    Area1 = abs(area1)+abs(area2);
    
    Edge = iner_edge(J,:);
    appose_node = zeros(Ne,2);
    for i = 1:Ne
        appose_node(i,1)  = Elem1(i,~((Elem1(i,:) == Edge(i,1)) + (Elem1(i,:) == Edge(i,2))));
        appose_node(i,2)  = Elem2(i,~((Elem2(i,:) == Edge(i,1)) + (Elem2(i,:) == Edge(i,2))));
    end
    
    Elem3 = [appose_node,Edge(:,1)];
    ve2 = node(Elem3(:,1),:)-node(Elem3(:,3),:);
    ve3 = node(Elem3(:,2),:)-node(Elem3(:,1),:);
    area3 = 0.5*(-ve3(:,1).*ve2(:,2) + ve3(:,2).*ve2(:,1));
    
    
     Elem4 = [appose_node,Edge(:,2)];
    ve2 = node(Elem4(:,1),:)-node(Elem4(:,3),:);
    ve3 = node(Elem4(:,2),:)-node(Elem4(:,1),:);
    area4 = 0.5*(-ve3(:,1).*ve2(:,2) + ve3(:,2).*ve2(:,1));
    
    Area2 = sort(abs([area3,area4]),2);
    Area2 = Area2(:,2);
    
    idx = Area1>(Area2 + eps);
    s(idx) = s(idx)+1;
    
           

end

