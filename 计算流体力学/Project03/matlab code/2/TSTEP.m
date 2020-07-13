function  Dt = TSTEP(rho,u,v,mu,DX,DY,a,K)
global gama Pr;
 a = a(2:end-1,2:end-1);
 rho = rho(2:end-1,2:end-1);
 u = u(2:end-1,2:end-1);
 v = v(2:end-1,2:end-1); 
 mu = mu(2:end-1,2:end-1); 
 v_1 = max(4*gama.*mu./(3*Pr.*rho));
 Dt= min(min(K.*(abs(u)/DX+abs(v)/DY+a.*sqrt(1/(DX^2)+1/(DY^2))...
                  +2*v_1*(1/(DX^2)+1/(DY^2))).^(-1)));
end