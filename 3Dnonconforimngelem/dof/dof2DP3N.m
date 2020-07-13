function [elem2Ddof3,face] = dof2DP3N(elem)
%%确定2D上p3-Nonconforming元的自由度,分别为面1，2，3和体上
%%为la1^2,la1*la2,la2^2;

 
NT = size(elem,1);
totalface = uint32([elem(:,[2 3]); elem(:,[3 1]); elem(:,[1 2])]);
totalfaces = sort(totalface,2);
[face,i2,j] = myunique(totalfaces);

%%
%面编号
l = length(i2);
faceL = 1:l;
faceL = [faceL;faceL+l;faceL+2*l];
faceL = faceL(:,j');

F=zeros(size(totalface));
Fs=zeros(size(totalface))+[1,2];
Fs = Fs';
for k=1:2
    F(:,k)=Fs((totalfaces==totalface(:,k))');
end

Fs = [F,F(:,1)+F(:,2)];
Fs = Fs';
Fss = zeros(size(totalface,1),3);

for k = 1:3
    Fss(:,k) =faceL(Fs==k);
end

elem2D3face = [Fss(1:NT,:),Fss(NT+1:2*NT,:),...
                        Fss(2*NT+1:3*NT,:)];
elem3ti = (1:NT)'+3*l;

elem2Ddof3=uint32([elem2D3face,elem3ti]);
end



