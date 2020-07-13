function [elem2dof,edge,bdDof] = dof2DP2N(elem)
%% DOFP2 dof structure for2D P2 Nonconforming finite element.
%


totalEdge = uint32(sort([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])],2));
[edge, i2, j] = myunique(totalEdge);
N = max(elem(:)); 
NT = size(elem,1);
NE = size(edge,1);
elem2edge = reshape(j,NT,3);
elemTi = (1:NT)'+(N+NE);
elem2dof = uint32([elem N+elem2edge elemTi]);
i1(j(3*NT:-1:1)) = 3*NT:-1:1; 
i1 = i1';
bdEdgeIdx = (i1 == i2);
isBdDof = false(N+2*NE,1);
isBdDof(edge(bdEdgeIdx,:)) = true;   % nodal 
idx = find(bdEdgeIdx);
isBdDof(N+idx) = true;
bdDof = find(isBdDof);