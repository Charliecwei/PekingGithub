function s = P2(Lam,X,R)
s = 0.5*sqrt(Lam)*(X(:,1).^2+X(:,2).^2-2*R^2);
end