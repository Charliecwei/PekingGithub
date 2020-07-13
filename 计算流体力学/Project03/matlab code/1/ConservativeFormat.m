function [x,U] = ConservativeFormat(pde,mesh)
%%Solve the equation u_t+f_x=0 in conservative form
global solver ;


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
    %max(abs(u))
   % maxvel
    t = t + mesh.h*mesh.nv/maxvel;
    lam = mesh.nv/maxvel;
    
    switch solver
        case 'LF'
            [F_l,F_r] = F_LFCr(Uo,pde.f,pde.a,lam,pde.gama);
%             F_l = F_LFCr(Uo(:,1:end-2),Uo(:,2:end-1),pde.f,lam,pde.gama);
%             F_r = F_LFCr(Uo(:,2:end-1),Uo(:,3:end),pde.f,lam,pde.gama);

    end
           
    Uo(:,2:end-1) = Uo(:,2:end-1) - lam*(F_r-F_l);
%    plot(xo,Uo(1,:));pause(0.1);
end

idx =(xo>=mesh.I(1))&(xo<=mesh.I(2));
U = Uo(:,idx);
x = xo(idx);

end


function [F_l,F_r] = F_LFCr(U,f,a,lam,gama)
%% the flux of LF scheme and correct
% Ur = U(:,2:end);
% Ul = U(:1:end-1);
% Fr = f(Ur,gama);
% Fl = f(Ul,gama);
F = f(U,gama);
Fl = F(:,1:end-1);
Fr = F(:,2:end);

Ul = U(:,1:end-1);
Ur = U(:,2:end);
Flux_LF = 0.5*(Fl+Fr)-0.5*(Ur-Ul)/lam;

Dletafplus = 0.5*(Fr-Fl) + 0.5*(Ur-Ul)/lam;
Deltafminus = 0.5*(Fr-Fl) - 0.5*(Ur-Ul)/lam;



A = a(0.5*(Ul+Ur),gama);

alphaDeltfplus = 0.25*(Dletafplus-lam*(A(:,:,1).*Dletafplus(1,:)...
                                 +A(:,:,2).*Dletafplus(2,:)+A(:,:,3).*Dletafplus(3,:)));
                             
alphaDeltfaminus = 0.25*(Deltafminus+lam*(A(:,:,1).*Deltafminus(1,:)...
                                 +A(:,:,2).*Deltafminus(2,:)+A(:,:,3).*Deltafminus(3,:)));


rplus = alphaDeltfplus(:,1:end-1)./alphaDeltfplus(:,2:end);
rminus = alphaDeltfaminus(:,2:end)./alphaDeltfaminus(:,1:end-1);




            
% rplus = [zeros(3,1)+1,rplus,zeros(3,1)+1];
% rminus  = [zeros(3,1)+1,rminus,zeros(3,1)+1];

rplus = [rplus,zeros(3,1)+1];
rminus = [zeros(3,1)+1,rminus];


% phisplus =  phi(rplus(:,1:end-1));
% phisminus = phi(rminus(:,2:end));

phisplus =  phi(rplus);
phisminus = phi(rminus);


FluxCor = Flux_LF + phisplus.*alphaDeltfplus - phisminus.*alphaDeltfaminus;
            


            

F_l =  FluxCor(:,1:end-1);
F_r =  FluxCor(:,2:end);

end



function  p = phi(r)
global phivar;

switch phivar
    case 'zero'
        p = zeros(size(r));
    case 'one'
      p = zeros(size(r))+1;
    case 'minmod'
        p = min(r,1);
        idx = r<0;
        p(idx) = 0;

    case 'Superbee'
           
            p = max(0,max(min(2*r,1),min(r,2)));

    case 'vanLeer'
        p = zeros(size(r));
        idx = r>0;
        p(idx) = 2*r(idx)./(1+r(idx));
end
end


