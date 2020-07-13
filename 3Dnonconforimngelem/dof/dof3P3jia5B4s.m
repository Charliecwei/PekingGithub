function [elem3dof3,face] = dof3P3jia5B4s(elem)
%%确定p3jia8B4的自由度,分别为面1，2，3，4和体上
%%为la1^2,la2^2,la3^2,la1la2,la1la3,la2la3;

 
NT = size(elem,1);
totalface = uint32([elem(:,[2 3 4]); elem(:,[3 4 1]); elem(:,[4 1 2]); ...
                    elem(:,[1 2 3 ])]);
totalfaces = sort(totalface,2);
[face,i2,j] = myunique(totalfaces);

%%
%面编号
l = length(i2);
faceL = 1:l;
faceL = [faceL;faceL+l;faceL+2*l;faceL+3*l;faceL+4*l;faceL+5*l];
faceL = faceL(:,j');

F=zeros(size(totalface));
Fs=zeros(size(totalface))+[1,2,3];
Fs = Fs';
for k=1:3
    F(:,k)=Fs((totalfaces==totalface(:,k))');
end

Fs = [F,F(:,1)+F(:,2)+1,F(:,1)+F(:,3)+1,F(:,2)+F(:,3)+1];
Fs = Fs';
Fss = zeros(size(totalface,1),6);

for k = 1:6
    Fss(:,k) =faceL(Fs==k);
end

elem3face = [Fss(1:NT,:),Fss(NT+1:2*NT,:),...
                        Fss(2*NT+1:3*NT,:),Fss(3*NT+1:4*NT,:)];
elem3ti = (1:NT)'+6*l;

elem3dof3=uint32([elem3face,elem3ti]);
end



