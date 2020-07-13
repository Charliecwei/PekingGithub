clear
global c_v R T_w mu_0 T_0 c_p Pr u_inf p_inf T_inf gama boundarycase;

i_max = 70;
j_max = 70;
MAXIT = 10000;
K = 0.6;
free_T = 1;


Ma_inf = 4.0;
LHORI = 0.00001;
a_inf = 340.28;
p_inf = 101325.0;
T_inf = 288.16;
 
T_w = T_inf*free_T;
gama = 1.4;
Pr = 0.71;
mu_0 = 1.789e-5;
T_0 = 288.16;
R = 287;



mu_inf = DYNVIS(T_inf);

V_inf = Ma_inf*a_inf;
u_inf = V_inf;
c_v = R/(gama-1);
c_p = gama*c_v;
rho_inf = p_inf/(R*T_inf);
Re_L = rho_inf*V_inf*LHORI/mu_inf;
e_inf = c_v*T_inf;

k_inf =  THERMAC(mu_inf);




delta = 5*LHORI/sqrt(Re_L);
LVERT = 5*delta;
Dx = LHORI/(i_max-1);
Dy = LVERT/(j_max-1);


%% 等温
boundarycase = 'isothermal';
u1 = zeros(i_max,j_max)+u_inf;
v1 = zeros(i_max,j_max);
p1 = zeros(i_max,j_max)+p_inf;
T1 = zeros(i_max,j_max)+T_inf;


[u1,v1,p1,T1] = Boundary(u1,v1,p1,T1);
rho1 = p1./(R.*T1);
e1 = c_v*T1;
mu1 = DYNVIS(T1);
k1 = THERMAC(mu1);

 %initial value
ITER1 = 1;
Is = 0;


while (ITER1<MAXIT)&&(Is==0)
    

    a = (gama.*R.*T1).^0.5;
    Dt = TSTEP(rho1,u1,v1,mu1,Dx,Dy,a,K); %Calculate the time step
    
   [rhos,u1,v1,p1,e1,T1,mu1,k1] = MAC(rho1,u1,v1,p1,e1,T1,mu1,k1,Dt,Dx,Dy); %Call the MAC algorithm
   
   Is = CONVER(rho1-rhos); %Check convergence
 
    rho1 = rhos;
  
    ITER1 = ITER1+1;
    
   
   
end



[deviation1,Out1] = MDOT(rho1,u1,Dy);%Check the subroutine for the conservation of mass






%% 绝热
boundarycase = 'adiabatic';
u2 = zeros(i_max,j_max)+u_inf;
v2 = zeros(i_max,j_max);
p2 = zeros(i_max,j_max)+p_inf;
T2 = zeros(i_max,j_max)+T_inf;


[u2,v2,p2,T2] = Boundary(u2,v2,p2,T2);
rho2 = p2./(R.*T2);
e2 = c_v*T2;
mu2 = DYNVIS(T2);
k2 = THERMAC(mu2);

 %initial value
ITER2 = 1;
Is = 0;


while (ITER2<MAXIT)&&(Is==0)
    

    a = (gama.*R.*T2).^0.5;
    Dt = TSTEP(rho2,u2,v2,mu2,Dx,Dy,a,K); %Calculate the time step
    
   [rhos,u2,v2,p2,e2,T2,mu2,k2] = MAC(rho2,u2,v2,p2,e2,T2,mu2,k2,Dt,Dx,Dy); %Call the MAC algorithm
   
   Is = CONVER(rho2-rhos); %Check convergence
 
    rho2 = rhos;
  %  sum(~isreal(u))
    ITER2 = ITER2+1;
   
   
end



[deviation2,Out2] = MDOT(rho2,u2,Dy);%Check the subroutine for the conservation of mass





















%OUTPUT;
%Generate a drawing data file
imax = i_max;
jmax = j_max;
y = 0:25/69:25;
figure
plot(p1(:,1)/p_inf,'-*');hold on
plot(p2(:,1)/p_inf,'-*');legend('等温壁','绝热壁')
ylim([0,6])
title('物面压力分布')
figure
plot(p1(imax,:)/p_inf,y,'-*')
hold on
plot(p2(imax,:)/p_inf,y,'-*')
legend('等温壁','绝热壁')
xlim([0,2.5])
title('边界层压力剖面')
figure
plot(T1(imax,:)/T_inf,y,'-*')
hold on
plot(T2(imax,:)/T_inf,y,'-*')
legend('等温壁','绝热壁')
xlim([0,4])
title('边界层内的温度剖面')


figure 
plot(u1(imax,:)/u_inf,y,'-*')
hold on
plot(u2(imax,:)/u_inf,y,'-*')
legend('等温壁','绝热壁')
title('边界层内的速度剖面')










