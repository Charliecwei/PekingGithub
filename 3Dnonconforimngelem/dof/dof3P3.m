function [elem3dof3,face,edge] = dof3P3(elem)
%%确定P3的自由度,分别为顶点1,2,3,4,面1,2;3,4和边1,2；1,3；
%%1,4;2,3;2,4;3,4

elemT = elem;

NT = size(elem,1);
totalface = uint32([elem(:,[2 3 4]); elem(:,[3 4 1]); elem(:,[4 1 2]); ...
                    elem(:,[1 2 3 ])]);
totalface = sort(totalface,2);
[face,~,j] = myunique(totalface);

elemF = reshape(j,NT,4)+max(abs(elem(:)));


totaledge = uint32([elem(:,[1,2]);elem(:,[1,3]);elem(:,[1,4]);...
                                         elem(:,[2,3]);elem(:,[2,4]);elem(:,[3,4])]);
totaledges = sort(totaledge,2);

[edge,~,j] = myunique(totaledges);

l =size(edge,1);

Eg = 1:l;
Eg = [Eg;Eg+l];
Eg = Eg(:,j');
eg = zeros(size(totaledge));
for jj = 1:2
    
    eg(:,jj) = Eg((totaledges==totaledge(:,jj))');
end

elemE = [eg(1:NT,:),eg(NT+1:2*NT,:),eg(2*NT+1:3*NT,:),eg(3*NT+1:4*NT,:),eg(4*NT+1:5*NT,:),eg(5*NT+1:6*NT,:)]+max(elemF(:));

elem3dof3 = [elemT,elemF,elemE];
