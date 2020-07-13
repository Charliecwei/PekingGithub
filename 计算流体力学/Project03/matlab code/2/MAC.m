function [rho,u,v,p,e,T,mu,k] = MAC(rho,u,v,p,e,T,mu,k,Dt,Dx,Dy)
% MAC scheme
global c_v R;


%% Build U
U_1 = rho;
U_2 = rho.*u;
U_3 = rho.*v;
E_t = rho.*(e+0.5*(power(u,2)+power(v,2)));
U_5 = E_t;

% forcast
state = 'forecast';
%% Build E
var = 'E3andE5';
tau_xx = TAUXX(mu,u,v,Dx,Dy,state);
tau_xy = TAUXY(mu,u,v,Dx,Dy,state,var);
q_x = QX(k,T,Dx,state);

E_1 = rho.*u;
E_2 = rho.*u.^2+p-tau_xx;
E_3 = rho.*u.*v-tau_xy;
E_5 = (E_t+p).*u-u.*tau_xx-v.*tau_xy+q_x;

%% Build F
var = 'F2andF5';
tau_yy = TAUYY(mu,u,v,Dx,Dy,state);
tau_xy =  TAUXY(mu,u,v,Dx,Dy,state,var);
q_y= QY(k,T,Dy,state);

F_1 = rho.*v;
F_2 = rho.*u.*v-tau_xy;
F_3 = rho.*v.^2+p-tau_yy;
F_5 = (E_t+p).*v-u.*tau_xy-v.*tau_yy+q_y;

%% forecast -- the forward difference
U_f_1 = U_1;
U_f_2 = U_2;
U_f_3 = U_3;
U_f_5 = U_5;

U_f_1(2:end-1,2:end-1) = U_1(2:end-1,2:end-1)- Dt/Dx*(E_1(3:end,2:end-1)-E_1(2:end-1,2:end-1))-Dt/Dy*(F_1(2:end-1,3:end)-F_1(2:end-1,2:end-1)); 
U_f_2(2:end-1,2:end-1) = U_2(2:end-1,2:end-1)- Dt/Dx*(E_2(3:end,2:end-1)-E_2(2:end-1,2:end-1))-Dt/Dy*(F_2(2:end-1,3:end)-F_2(2:end-1,2:end-1)); 
U_f_3(2:end-1,2:end-1) = U_3(2:end-1,2:end-1)- Dt/Dx*(E_3(3:end,2:end-1)-E_3(2:end-1,2:end-1))-Dt/Dy*(F_3(2:end-1,3:end)-F_3(2:end-1,2:end-1)); 
U_f_5(2:end-1,2:end-1) = U_5(2:end-1,2:end-1)- Dt/Dx*(E_5(3:end,2:end-1)-E_5(2:end-1,2:end-1))-Dt/Dy*(F_5(2:end-1,3:end)-F_5(2:end-1,2:end-1)); 


[u,v,p,T] = DECODE(U_f_1,U_f_2,U_f_3,U_f_5); %%git rho,u,v,p of inner point



%% boundary condition

[u,v,p,T] = Boundary(u,v,p,T);

rho = p./(R.*T);
e = c_v*T;
mu = DYNVIS(T);
k = THERMAC(mu);

state = 'correct';
%% Build E
var = 'E3andE5';
tau_xx = TAUXX(mu,u,v,Dx,Dy,state);
tau_xy = TAUXY(mu,u,v,Dx,Dy,state,var);
q_x = QX(k,T,Dx,state);

E_t = rho.*(e+0.5*(power(u,2)+power(v,2)));
E_1 = rho.*u;
E_2 = rho.*u.^2+p-tau_xx;
E_3 = rho.*u.*v-tau_xy;
E_5 = (E_t+p).*u-u.*tau_xx-v.*tau_xy+q_x;

%% Build F
var = 'F2andF5';
tau_yy = TAUYY(mu,u,v,Dx,Dy,state);
tau_xy =  TAUXY(mu,u,v,Dx,Dy,state,var);
q_y= QY(k,T,Dy,state);


F_1 = rho.*v;
F_2 = rho.*u.*v-tau_xy;
F_3 = rho.*v.^2+p-tau_yy;
F_5 = (E_t+p).*v-u.*tau_xy-v.*tau_yy+q_y;

%% Correct -- backward difference
U_1(2:end-1,2:end-1) = 0.5*((U_1(2:end-1,2:end-1)+U_f_1(2:end-1,2:end-1))...
             - Dt/Dx*(E_1(2:end-1,2:end-1)-E_1(1:end-2,2:end-1))...
              -Dt/Dy*(F_1(2:end-1,2:end-1)-F_1(2:end-1,1:end-2)));

U_2(2:end-1,2:end-1) = 0.5*((U_2(2:end-1,2:end-1)+U_f_2(2:end-1,2:end-1))...
             - Dt/Dx*(E_2(2:end-1,2:end-1)-E_2(1:end-2,2:end-1))...
              -Dt/Dy*(F_2(2:end-1,2:end-1)-F_2(2:end-1,1:end-2)));

U_3(2:end-1,2:end-1) = 0.5*((U_3(2:end-1,2:end-1)+U_f_3(2:end-1,2:end-1))...
             - Dt/Dx*(E_3(2:end-1,2:end-1)-E_3(1:end-2,2:end-1))...
              -Dt/Dy*(F_3(2:end-1,2:end-1)-F_3(2:end-1,1:end-2)));

U_5(2:end-1,2:end-1) = 0.5*((U_5(2:end-1,2:end-1)+U_f_5(2:end-1,2:end-1))...
             - Dt/Dx*(E_5(2:end-1,2:end-1)-E_5(1:end-2,2:end-1))...
              -Dt/Dy*(F_5(2:end-1,2:end-1)-F_5(2:end-1,1:end-2)));



[u,v,p,T] = DECODE(U_1,U_2,U_3,U_5); %%git rho,u,v,p;

%% boundary condition

[u,v,p,T] = Boundary(u,v,p,T);

rho = p./(R.*T);
e = c_v*T;
mu = DYNVIS(T);
k = THERMAC(mu);



 
end




