function [x,U] = ConservativeFormat(pde,option,mesh)
%%Solve the equation u_t+f_x=0 in conservative form



switch pde.bdtype
    case 'Inf'
        xo = mesh.I(1):mesh.h:mesh.I(2);
        Uo = pde.uo(xo);
         rho =( Uo(1,1:end-1)+Uo(1,2:end))/2;
        u = (Uo(2,1:end-1)+Uo(2,2:end))./(2*rho);
        E = (Uo(3,1:end-1)+Uo(3,2:end))/2;
        p = (pde.gama-1)*(E-0.5*rho.*u.^2);
         maxvel = max(abs(u)+sqrt(pde.gama*p./rho));
         Tall = maxvel*mesh.T+1;
        Io = [mesh.I(1)-Tall,mesh.I(2)+Tall];
%           Io = [-10,10];
    case 'Fixed'
        Io = mesh.I;
end

xo = Io(1):mesh.h:Io(2);
Uo = pde.uo(xo);
t = 0;
while t<mesh.T
    rho =( Uo(1,1:end-1)+Uo(1,2:end))/2;
    u = (Uo(2,1:end-1)+Uo(2,2:end))./(2*rho);
    E = (Uo(3,1:end-1)+Uo(3,2:end))/2;
    p = (pde.gama-1)*(E-0.5*rho.*u.^2);

%    maxvel = max(abs(Uo(2,:)./Uo(1,:))+sqrt(pde.gama*p./Uo(1,:)));
% 
%     rho =Uo(1,:);
%     u = Uo(2,:)./rho;
%     E = Uo(3,:);
%     p = (pde.gama-1)*(E-0.5*rho.*u.^2);
% 
% 
     maxvel = max(abs(u)+sqrt(pde.gama*p./rho));
   % maxvel
    t = t + mesh.h*mesh.nv/maxvel;
    lam = mesh.nv/maxvel;
    
    switch option.solver
        case 'LF'
            F_l = F_LF(Uo(:,1:end-2),Uo(:,2:end-1),pde.f,lam,pde.gama);
            F_r = F_LF(Uo(:,2:end-1),Uo(:,3:end),pde.f,lam,pde.gama);
        case 'LW'
            F_l = F_LW(Uo(:,1:end-2),Uo(:,2:end-1),pde.f,pde.a,lam,pde.gama);
            F_r = F_LW(Uo(:,2:end-1),Uo(:,3:end),pde.f,pde.a,lam,pde.gama);
%               u = Uo(:,1:end-2)';
%               v = Uo(:,2:end-1)';
%               F_l =  EulerLW(u,v,pde.gama,lam);
%               u = Uo(:,2:end-1)';
%               v = Uo(:,3:end)';
%                F_r =  EulerLW(u,v,pde.gama,lam);
%                F_l = F_l';
%                F_r = F_r';

    end
    Uo(:,2:end-1) = Uo(:,2:end-1) - lam*(F_r-F_l);
 %  plot(xo,Uo(1,:));pause(0.1)
end

idx =(xo>=mesh.I(1))&(xo<=mesh.I(2));
U = Uo(:,idx);
x = xo(idx);

end


function Flux = F_LF(U,V,f,lam,gama)
%% the flux of LF scheme
Flux = 0.5*(f(U,gama)+f(V,gama))-0.5*(V-U)/lam;
end


function Flux = F_LW(U,V,f,a,lam,gama)
%% the flux of LW scheme
A = a((U+V)/2,gama);
Fv = f(V,gama);
Fu = f(U,gama);
F = Fv - Fu;
Flux = 0.5*(Fu+Fv)-0.5*lam*(A(:,:,1).*F(1,:)+A(:,:,2).*F(2,:)+A(:,:,3).*F(3,:));
end

                    
% function [numflux] = EulerLW(u,v,gamma,lambda)
% % function [numflux] = EulerLW(u,v,gamma,lambda,maxvel);
% % Purpose : Evaluate Lax Wendroff numerical flux for the Euler equations
% % Compute flux for u
% r = u(:,1); ru = u(:,2); E = u(:,3); pu = (gamma-1)*(E - 0.5*ru.^2./r);
% fu = [ru (ru.^2./r+pu) (E+pu).*ru./r];
% % Compute flux for v
% r = v(:,1); ru = v(:,2); E = v(:,3); pv = (gamma-1)*(E - 0.5*ru.^2./r); 
% fv = [ru (ru.^2./r+pv) (E+pv).*ru./r];
% % Evaluate numerical flux
% fw = fv-fu; w= (u+v)/2; rw =w(:,1); ruw =w(:,2);
% Ew=w(:,3); uw=ruw./rw; wL= length(rw);
% A21 = -(3-gamma)/2*uw.^2; A22 = (3-gamma)*uw; A23 = (gamma-1)*ones(wL,1) ; A31 = -gamma*Ew.*uw./rw + (gamma-1)*uw.^3;
% A32 = gamma*Ew./rw - 3*(gamma-1)/2*uw.^2; A33 = gamma*uw;
% LFvec = zeros(wL,3); LFvec(:,1) = fw(:,2);
% LFvec(:,2) = A21.*fw(:,1)+A22.*fw(:,2)+A23.*fw(:,3); LFvec(:,3) = A31.*fw(:,1)+A32.*fw(:,2)+A33.*fw(:,3); numflux = ((fu+fv) - lambda*LFvec)/2;
% return            
% end