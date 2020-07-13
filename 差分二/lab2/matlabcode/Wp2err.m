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
% idx = ~isnan(s);
 Wp2 = power(sum(s.*w),2/P);
 
end

 
 