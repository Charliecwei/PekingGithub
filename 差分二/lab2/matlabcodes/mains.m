
%% Initialize
global u0 n node N0 ;
pde = example1;
N = [2,4,8,16,32];
P = [1,2]; 
R = 1.5;
Lam = 25;
Tol = 1e-5;
Errinf = zeros(length(N),1);
Wp2Err = zeros(length(N),length(P));

Errinf_Rate = zeros(length(N),1);
Wp2Err_Rate = zeros(length(N),length(P));
InteK = zeros(length(N),1);
Hises = cell(length(N),1);

%% Compution

for k = 1:length(N)
    fprintf('当前 n = %d, 计算进度\n',N(k));
    fprintf('00%%--20%%--40%%--60%%--80%%-100%%\n')
    n = N(k);
     [InteK(k),Hises{k}] = Oliker_Prussner(pde,R,Lam,Tol); 
    
       % subtitle('h=%3.f',2/N(k));

     %% exact solution
    exactu = pde.u(node);
    exactu = exactu(1:N0);
    %% err 
    err = exactu - u0;
  
    % Inf err
    Errinf(k) = norm(abs(err),inf);
    if k>1
        Errinf_Rate(k) = log( Errinf(k-1)/ Errinf(k))/log(N(k)/N(k-1));
    end
    
    % Wp2 err
    for ep = 1:length(P)
           Wp2Err(k,ep) = Wp2err(err,P(ep));
           
           if k>1
                Wp2Err_Rate(k,ep) = log(Wp2Err(k-1,ep)/Wp2Err(k,ep))/log(N(k)/N(k-1));
          end
    end
    
end
    





 %% 打印结果
 Print(N,P,Errinf,Errinf_Rate,Wp2Err,Wp2Err_Rate,InteK,Hises);


 

 
 
function  [InteK,Hises] = Oliker_Prussner(pde,R,Lam,Tol)

global node elem iner_edge  node_iner_edge_idx node_elem_idx iner_edge_elem_idx ;
global u0 N0 ubd  n;

nx = n;
ny = n;
hx = 2/nx;
hy = 2/ny;
%h = max(hx,hy);

x = (-1:hx:1)';
y = (-1:hy:1)';
node = [kron(ones(ny+1,1),x),kron(y,ones(nx+1,1))];
idx = ~or(or(node(:,1)==1,node(:,1)==-1),or(node(:,2)==1,node(:,2)==-1));
node = [node(idx,:);node(~idx,:)];
N0 = sum(idx);
%node_iner = node(1:N0,:);


Ndof = size(node,1);

elem = uint16(sort(delaunay(node(:,1),node(:,2)),2));


totalEdge = uint16([elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])]);
totalEdges = sort(totalEdge,2);
[edge, i2, J] = unique(totalEdges,'rows','legacy');


NT = size(elem,1);
i1(J(3*NT:-1:1)) = 3*NT:-1:1; 
i1 = i1';
IntEdgeIdx = find(~(i1 == i2));
iner_edge = uint16(edge(IntEdgeIdx,:));
NE =  length(IntEdgeIdx);
iner_edge_elem_idx = uint16(zeros(NE,2));
J =  reshape(J,[],3);
for i = 1:length(IntEdgeIdx)
    idx = find(sum(J == IntEdgeIdx(i),2)==1)';
    iner_edge_elem_idx(i,:) = uint16(idx);
end

node_iner_edge_idx = cell(Ndof,1);
node_elem_idx = cell(Ndof,1);
for i = 1:Ndof
    node_iner_edge_idx{i} = uint16(find(sum(iner_edge==i,2)))';
    node_elem_idx{i} = uint16(find(sum(elem==i,2)))';
end

for i = 1:N0
    I = i;
    Elem = elem(node_elem_idx{I},:)';
    edge_ap = sort(reshape(Elem(Elem~=I),2,[]),1)';
    ev = (node(edge_ap(:,1),:)+node(edge_ap(:,2),:))-2*node(I,:);
    k = convhull(ev(:,1),ev(:,2));
    node_elem_idx{I} = node_elem_idx{I}(k(1:end-1));
end

%%




%% Compute f

 [~,areas] = gradbasis(node,elem);
NT = size(elem,1);
[lambda,weight] = quadpts(4);

phi = lambda;                 % linear bases
nQuad = size(lambda,1);
bt = zeros(NT,3);
for p = 1:nQuad
		% quadrature points in the x-y-z coordinate
		pxyz = lambda(p,1)*node(elem(:,1),:) ...
			 + lambda(p,2)*node(elem(:,2),:) ...
			 + lambda(p,3)*node(elem(:,3),:);
		fp = pde.f(pxyz);
        for js = 1:3
            bt(:,js) = bt(:,js) + weight(p)*phi(p,js)*fp;
        end
end
   
    bt = bt.*repmat(areas,1,3);
    f = accumarray(elem(:),bt(:),[Ndof 1]);
    f = f(1:N0);
    
    if min(f)<0
         fprintf('Wrong! f have minus\n');
         close
    end
        
    clear pxyz bt
    
%%    



%% initial u0


u0 = P2(Lam,node,R);




%% let u0(bd) = ubd and inucled mesh
  ubd = pde.g_D(node(N0+1:end,:));

if ~isempty(find(uint16(u0(N0+1:end) < ubd)==0, 1))
    fprintf('Wrong! Lam or R small\n');
     close
end



for i = 1:length(ubd)
   u0(i+N0) = ubd(i);
   Inducled_mesh(i+N0);
end


%% check the convex 

jump = Jump(1:NE);
idx = find(jump<0, 1);
if ~isempty(idx)
    fprintf('Wrong! u0 at the mesh is not convex \n');
    close
end



%% Compute the initial S  
S = zeros(N0,1);
for i = 1:N0
    S(i) =  Conhulls(i);
%     s =  Conhull(i);
%     if s~=S(i)
%         stop
%     end
    if S(i)< f(i)
         fprintf('Wrong! Lam  small\n');
         close
    end
end

%%  Iterative update u0

es = max(S-f); % es is the max for current all node S-f

Es = log(es/Tol);
KT = 1;
InteK = 0;
Hises = [];
while es>Tol
    InteK = InteK+1;
    Hises = [Hises;es];
    for i = 1:N0
        I = i;
         S(i) = Conhulls(i);
         e = S(i) - f(i);       % e is the current update node I S-f
         while e>1e-9

                thresold =Thresold(I);
                [a,b] = min(thresold);
                
                au = a(1);
                
                u0(I) = u0(I)+au;
                S(i) = Conhulls(I);
                 e = S(i) - f(i);

                 if e<0
                     u0(I) = u0(I) - au;
                      delta = Find_deltas(I,au,f(i));
                     u0(I) = u0(I)+delta;
                      S(I) = Conhull(I);
                        e = S(I) - f(I);
                 else
                     Flips(node_iner_edge_idx{I}(b(1)));
                 end
         end
    end
                     
                
    
    for i = 1:N0

        S(i) =  Conhulls(i);
%          s =  Conhull(i);
%     if s~=S(i)
%         stop
%     end
    if S(i)< f(i)-1e-10
           fprintf('Wrong! The update has bug\n');
         close
    end
    end
    es = max(S-f);
    if (log(es/Tol)<KT*Es)
        KT = KT-1/25;
        fprintf('.')
    end
end

fprintf('\n')  


u0 = u0(1:N0);      

end









function Wp2 = Wp2err(err,P)
global node elem N0 
 %% compute Wp2 err
 e = [1,0;0,1;1,1;1,-1];
 s = zeros(N0,1);
 Ndof = size(node,1);
 for k = 1:size(e,1)
     s = s + abs(Delta_e(err,e(k,:))).^P;
 end
 
 [~,areas] = gradbasis(node,elem);
 bt = repmat(areas,1,3);
 w = accumarray(elem(:),bt(:),[Ndof 1]);
 w = w(1:N0);
 Wp2 = power(sum(s.*w),1/P);
 
end

 
 