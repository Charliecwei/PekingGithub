function s =LocaltoGlobalidx(i,j,N)
%%The position of the local interior point (i,j) in the global interior point
%is pulled up in the y direction of size N 
%(that means there are N-1 inner point on each x or y diection)

s = (j-1).*(N-1)+i;
end


