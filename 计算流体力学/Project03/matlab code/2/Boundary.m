function [u,v,p,T] = Boundary(u,v,p,T)
%% git the boundary condition
global boundarycase;
global  T_w  u_inf p_inf T_inf;
% The first kind of boundary point--leading edge points
u(1,1) = 0;
v(1,1) = 0;
p(1,1) = p_inf;
T(1,1) = T_inf;

%The second kind of boundary point--inflow boundary and upper boundary
u(1,2:end) = zeros(size(u(1,2:end)))+u_inf;
v(1,2:end) = zeros(size(v(1,2:end)))+0;
p(1,2:end) = zeros(size(p(1,2:end)))+p_inf;
T(1,2:end) = zeros(size(T(1,2:end)))+T_inf;

u(2:end,end) = zeros(size(u(2:end,end))) + u_inf;
v(2:end,end) = zeros(size(v(2:end,end))) + 0;
p(2:end,end) = zeros(size(p(2:end,end)))+p_inf;
T(2:end,end) = zeros(size(T(2:end,end)))+T_inf;

%The third kind of boundary point --the surface

u(2:end,1) = zeros(size(u(2:end,1))) + 0;
v(2:end,1) = zeros(size(v(2:end,1))) + 0;
p(2:end,1) = 2*p(2:end,2) - p(2:end,3);
switch boundarycase
    case 'isothermal'
        T(2:end,1) = zeros(size(T(2:end,1))) + T_w;
    case 'adiabatic'
        T(2:end,1) = T(2:end,2);
end

%The fourth kind of boundary point --the flow boundary
u(end,2:end-1) = 2*u(end-1,2:end-1) - u(end-2,2:end-1);
v(end,2:end-1) = 2*v(end-1,2:end-1) - v(end-2,2:end-1);
p(end,2:end-1) = 2*p(end-1,2:end-1) - p(end-2,2:end-1);
T(end,2:end-1) = 2*T(end-1,2:end-1) - T(end-2,2:end-1);

end
