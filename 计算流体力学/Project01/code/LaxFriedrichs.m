function [U,CourantNum] = LaxFriedrichs(Uo,h,tau,a,T,bdtype)
%   the space step length is h, 
%   the time step length is tau, 
%   the termination time is T, 
%   the initial function data is Uo, 
%    hyperbolic equation u_t+au_x=0；
%   We assume a is a constant or a function about u.
%   bdtype is the boundary type, for Periodic boundary or Fixed boundary
%%
lam = tau/h;
N = T/tau;
    switch bdtype
        case 'Periodic'
            if isreal(a)
                 CourantNum = 0;
                 for n = 1:N
                    U = OnestepLaxFriedrichsPer(Uo,lam,a);
                    Uo = U;
                    CourantNum = max([CourantNum abs(a)*lam]);
                 end
            else
                  f = a.f;
                  a = a.a;
                 CourantNum = 0;
                 for n = 1:N
                    Flux = f(Uo); %通量
                    as = a(0.5*(Uo(1:end-1)+Uo(2:end)));
                    U = OnestepConLFPer(Uo,Flux,lam); %守恒型
                    Uo = U;
                    CourantNum = max([CourantNum abs(as)*lam]);
%                   ylim([-0.5,2]);
%                    plot(0:h:2*pi,U); pause(0.1); 
                 end
            end

        case 'Fixed'
                if isreal(a)
                     CourantNum = 0;
                     for n = 1:N
                        U = OnestepLaxFriedrichsFix(Uo,lam,a);
                        Uo = U;
                        CourantNum = max([CourantNum abs(a)*lam]);
                     end
                else
                      f = a.f;
                      a = a.a;
                     CourantNum = 0;
                     for n = 1:N
                        Flux = f(Uo); %通量
                        as = a(0.5*(Uo(1:end-1)+Uo(2:end)));
                        U = OnestepConLFFix(Uo,Flux,lam); %守恒型
                        Uo = U;
                        CourantNum = max([CourantNum abs(as)*lam]);
    %                   ylim([-0.5,2]);
    %                    plot(0:h:2*pi,U); pause(0.1); 
                     end
                 end
    end
end






function U = OnestepLaxFriedrichsPer(uo,lam,a)
%% One step for  CIR scheme
%uo is the discrete grid point value
%a is the discrete coefficient of the corresponding point or constant
%  Periodic boundary
%%
        L = length(uo);
        U = zeros(1,L);
            
            U(1) = 0.5*(uo(2)+uo(L-1)) - 0.5*a*lam*(uo(2)-uo(L-1));
            U (2:L-1) = 0.5*(uo(3:L) + uo(1:L-2)) - 0.5*a*lam*(uo(3:L)-uo(1:L-2));
            U(L) = 0.5*(uo(2)+uo(L-1)) - 0.5*a*lam*(uo(2)-uo(L-1));

end





function U = OnestepLaxFriedrichsFix(uo,lam,a)
%% One step for  CIR scheme
%uo is the discrete grid point value
%a is the discrete coefficient of the corresponding point or constant
%  Periodic boundary
%%
        L = length(uo);
        U = zeros(1,L);
        U(1) = uo(1);      
        U (2:L-1) = 0.5*(uo(3:L) + uo(1:L-2)) - 0.5*a*lam*(uo(3:L)-uo(1:L-2));
%       for i = 2:L-1
%           U(i) = 0.5*(uo(i+1)+uo(i-1)) - 0.5*a*lam*(uo(i+1)-uo(i-1));
%       end
        U(L) = uo(L);

end
    
    
function U = OnestepConLFPer(uo,Flux,lam)
    L = length(uo);
    U = zeros(1,L);
    U(1) = 0.5*(uo(2)+uo(L-1)) - 0.5*lam*(Flux(2)-Flux(L-1));
    U(2:L-1) = 0.5*(uo(3:L)+uo(1:L-2)) - 0.5*lam*(Flux(3:L) - Flux(1:L-2));
    U(L) = U(1);
end


   
function U = OnestepConLFFix(uo,Flux,lam)
    L = length(uo);
    U = zeros(1,L);
    U(1) = uo(1);
    U(2:L-1) = 0.5*(uo(3:L)+uo(1:L-2)) - 0.5*lam*(Flux(3:L) - Flux(1:L-2));
    U(L) = uo(L);
end
    
    
