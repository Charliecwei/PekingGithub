function [i,j] =GlobaltoLocalidx(s,N)
%%The position of the  in the global interior point to local interior point (i,j)
%is pulled up in the y direction of size N 
%(that means there are N-1 inner point on each x or y diection)

j = ceil(s/(N-1));
i = s-(j-1)*(N-1);
end