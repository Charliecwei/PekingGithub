function [elem3dof2,face,bdface,IdpN,elemface,Fo] = dof3P2jia(elem)
%%ȷ��p2jia�����ɶ�,�ֱ�Ϊ��1��2��3��4������


 
NT = size(elem,1);
totalface = uint32([elem(:,[2 3 4]); elem(:,[3 4 1]); elem(:,[4 1 2]); ...
                    elem(:,[1 2 3 ])]);
totalfaces = sort(totalface,2);
[face,i2,j] = myunique(totalfaces);

%%
%����
l=length(i2);
faceL=1:l;
faceL=[faceL;faceL+l;faceL+2*l];
faceL=faceL(:,j');
%faceLL=faceL(:);
F=zeros(size(totalface));
%F1=zeros(size(totalface));
for k=1:3
    F(:,k)=faceL((totalfaces==totalface(:,k))');
  %  g=(totalfaces==totalface(:,k))';
   % g=g(:);
  %  F1(:,k)=faceLL(g);
   % sum(F1(:,k)~=F(:,k))
end


elem3face = [F(1:NT,:),F(NT+1:2*NT,:),F(2*NT+1:3*NT,:),F(3*NT+1:4*NT,:)];
elem3ti=1:NT;
elem3ti=elem3ti'+3*l;

elem3dof2=uint32([elem3face,elem3ti]);













if nargout >=3
    %%����߽��� bdace(Idp(j)+1:Idp(j+1),:)��ʾ��j���棨�Ե�Ԫ��˵�ģ� Ϊ�߽���
    Nfall = length(totalface);
    i1(j(Nfall:-1:1)) = Nfall:-1:1;
    i1 = i1';
    bdFlag = zeros(Nfall,1,'uint8');
    bdFaceidx = i1(i1==i2);   
   bdFlag(bdFaceidx) = true;   
   bdFlag = reshape(bdFlag,NT,4);
   
  
   idx = 0;
   idp = 0;
   IdpN =zeros(5,1);
   for jj = 1:4
       Totalface = totalfaces(idx+1:idx+NT,:);
       g = Totalface(bdFlag(:,jj)==1,:);
       
       l=size(g,1);
       if l ~= 0
                 bdface(idp+1:idp+l,:) = g;
       end
       idx = idx +NT;
       idp = idp +l;
       IdpN(jj+1) = idp;
   end

elemface = reshape(j,NT,4);



faceL = zeros(3,size(totalface,1));
faceL(1,:) = 1;
faceL(2,:) = 2;
faceL(3,:) = 3;
Fo = zeros(size(totalface));

for k=1:3
    Fo(:,k)=faceL((totalface==totalfaces(:,k))');
end





end
