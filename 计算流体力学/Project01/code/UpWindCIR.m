function [U, CourantNum]= UpWindCIR(Uo,h,tau,a,T,bdtype)
%   the space step length is h, 
%   the time step length is tau, 
%   the termination time is T, 
%   the initial function data is Uo, 
%    hyperbolic equation u_t+au_x=0(a is a constant) or u_t+f(u)_x=0 (f is
%    a function of u and a = f'(u));
%   bdtype is the boundary type, for Periodic boundary or Fixed boundary
%%
lam = tau/h;
N = T/tau;
    switch bdtype
        case 'Periodic'
            if isreal(a)
                for n = 1:N
                    U = OnestepCIRPer(Uo,lam,a);
                    Uo = U;
                end
                  CourantNum = abs(a)*lam;
            else %% a is a struct
                f = a.f;
                a = a.a;
                 CourantNum = 0;
                 for n = 1:N
                    Flux = f(Uo); %通量
                    as = a(0.5*(Uo(1:end-1)+Uo(2:end)));
%                    as = zeros(1,length(Uo)-1);
%                    umid = Uo(2:end) - Uo(1:end-1);
%                    idx = (umid==0);
%                    as(idx) = a(Uo(idx));
%                    idx = find(umid);
%                    as(idx) = (f(Uo(idx+1))-f(Uo(idx)))/Uo(idx);
                    U = OnestepConCIRPer(Uo,Flux,lam,as); %守恒型
                    Uo = U;
                    CourantNum = max([CourantNum abs(as)*lam]);
                 end
            end

        case 'Fixed'
            if isreal(a)
                for n = 1:N
                    U = OnestepCIRFix(Uo,lam,a);
                    Uo = U;
      %               plot(-2:h:4,U);pause(0.1);
                end
               
                  CourantNum = abs(a)*lam;
            else
                 f = a.f;
                 a = a.a;
                 CourantNum = 0;
                 for n = 1:N
                    Flux = f(Uo); %通量
                    as = a(0.5*(Uo(1:end-1)+Uo(2:end)));
                    U = OnestepConCIRFix(Uo,Flux,lam,as); %守恒型
                    Uo = U;
                    CourantNum = max([CourantNum abs(as)*lam]);
            %         plot(-25.5:h:25.5,U); pause(0.1); 
                 end
            end
    end
end






function U = OnestepCIRPer(uo,lam,a)
%% One step for  CIR scheme
%uo is the discrete grid point value
%a is the discrete coefficient of the corresponding point or constant
%  Periodic boundary
%%
        L = length(uo);
        U = zeros(1,L);
        if length(a)==1
            if a>0
                U(2:L) = uo(2:L)-a*lam*(uo(2:L)-uo(1:L-1));
                U(1) = U(L);
            elseif a<0
                U(1:L-1) = uo(1:L-1)-a*lam*(uo(2:L)-uo(1:L-1));
                U(L) = U(1);
            else
                U(1:L) = uo(1:L);
            end
        else
            if a(1)>0
                U(1) = uo(1) - a(1)*lam*(uo(1)-uo(L-1));
            elseif a(1)<0
                U(1) = uo(1) - a(1)*lam*(uo(2)-uo(1));
            else
                U(1) = uo(1);
            end
            as = a(2:end-1);

            idx=find((as>0))+1;
            U(idx) = uo(idx)-as(idx)*lam*(uo(idx)-uo(idx-1));

            idx = find((as<0))+1;
            U(idx) = uo(idx)-as(idx)*lam*(uo(idx+1)-uo(idx));

            idx = find((as==0))+1;
            U(idx) = uo(idx);
            if a(L)>0
                U(L) = uo(L) - a(L)*lam*(uo(L)-uo(L-1));
            elseif a(L)<0
                U(L) = uo(L) - a(L)*lam*(uo(2)-uo(L));
            end
        end
end





function U = OnestepCIRFix(uo,lam,a)
%% One step for  CIR scheme
%uo is the discrete grid point value
%a is the discrete coefficient of the corresponding point or constant
%  Periodic boundary
%%
        L = length(uo);
        U = zeros(1,L);
        U(1) = uo(1);
            if a>0
%                for i = 2:L-1
%                    U(i) = uo(i) - a*lam*(uo(i)-uo(i-1));
%                end
                U(2:L-1) = uo(2:L-1)-a*lam*(uo(2:L-1)-uo(1:L-2));
            elseif a<0
                U(2:L-1) = uo(2:L-1)-a*lam*(uo(3:L)-uo(2:L-1));
            else
                U(2:L-1) = uo(2:L-1);
            end
        U(L) = uo(L);

end
    
    


function U = OnestepConCIRPer(uo,Flux,lam,a)
%% The conservative difference scheme solves the periodic boundary
a = abs(a);
L = length(uo);
U = zeros(1,L);
U(1) = uo(1) - 0.5*lam*(Flux(2) - Flux(L-1))...
                + 0.5*lam*a(1)*(uo(2)-uo(1))...
                -0.5*lam*a(L-1)*(uo(1)-uo(L-1));
           
U(2:L-1) = uo(2:L-1) - 0.5*lam*(Flux(3:L) - Flux(1:L-2)) ...
                      + 0.5*lam.*a(2:L-1).*(uo(3:L)-uo(2:L-1))...
                       -0.5*lam.*a(1:L-2).*(uo(2:L-1)-uo(1:L-2));

            
U(L) = uo(L) - 0.5*lam*(Flux(2)-Flux(L-1))...
                      +0.5*lam*a(1)*(uo(2)-uo(L))...
                      -0.5*lam*a(L-1)*(uo(L)-uo(L-1));

end

    
    

function U = OnestepConCIRFix(uo,Flux,lam,a)
%% The conservative difference scheme solves the periodic boundary
a = abs(a);
L = length(uo);
U = zeros(1,L);
U(1) = uo(1);
           
U(2:L-1) = uo(2:L-1) - 0.5*lam*(Flux(3:L) - Flux(1:L-2)) ...
                      + 0.5*lam.*a(2:L-1).*(uo(3:L)-uo(2:L-1))...
                       -0.5*lam.*a(1:L-2).*(uo(2:L-1)-uo(1:L-2));
                   
U(L) = uo(L);

end