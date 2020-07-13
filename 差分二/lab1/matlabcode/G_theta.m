function G = G_theta(N)
%%Find the direction of the size of N template
% G(i,1),G(i,2) are the lengths in the x and y directions for i-th theta
if N ==1
    G = [1,0;1,1];
else
    Gs = [];
    l = 0;
    %% 0<theta<pi/4
    for i = 2:N
        for j = 1:i-1
            l = l+1;
            k = gcd(i,j);
            Gs(l,:) = [i/k,j/k];
        end
    end
    Gs =  myunique(Gs);
    G = [[1,0];Gs;[1,1];[Gs(:,2),Gs(:,1)]];
end
