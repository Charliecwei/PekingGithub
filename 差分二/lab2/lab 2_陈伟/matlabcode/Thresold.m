function  thresoldH =Thresold(I)
%%find all  thresold value for point I and return it's min thresold,
%%adjacent node,appose node, Convex judge the edge can be change or not. 1
%%can be change 0 can't
global node elem iner_edge node_iner_edge_idx iner_edge_elem_idx u0 N0

node_edge_idx =  node_iner_edge_idx{I};
if isempty(node_edge_idx)
    thresoldH = inf;

else
   L_ne = length(node_edge_idx);
   thresoldH = zeros(L_ne,1);

    for i = 1:L_ne
        n_edgeI = node_edge_idx(i);
        Edge =  iner_edge(n_edgeI,:);
        ad_node = setdiff(Edge,I);
        Elem = elem(iner_edge_elem_idx(n_edgeI,:),:);
        ap_node = setdiff(Elem,Edge);
        X = [node([ad_node;ap_node],:),u0([ad_node;ap_node])]';
     %   size(X)
        A = X(1:2,2:3)-X(1:2,1);
        if rank(A) == 1
           thresoldH(i) = inf;
        else
            b = node(I,:)'-X(1:2,1);
            x = A\b;
            thresoldH(i) =X(3,1)+dot(X(3,2:3)-X(3,1),x) - u0(I);
            if and(thresoldH(i)<0,I<=N0)
                thresoldH(i) = inf;
            end
        end
            

            

    end
end

