function [deviation,Out] = MDOT(rho,u,Dy)
rho = rho';
u = u';
u_in = rho(:,1).*u(:,1);
u_out = rho(:,end).*u(:,end);

influx = sum(0.5*(u_in(2:end)+u_in(1:end-1)))*Dy;
outflux = sum(0.5*(u_out(2:end)+u_out(1:end-1)))*Dy;

net_inflow = influx - outflux;

deviation =  abs(net_inflow)/(influx+outflux);

if deviation <0.01
    Out = 'valid';
else
    Out = 'invalid';
    
end
end

