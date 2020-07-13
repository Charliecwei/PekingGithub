function mit18086_navierstokeshw(example)
%MIT18086_NAVIERSTOKES
%    Solves the incompressible Navier-Stokes equations in a
%    rectangular domain with prescribed velocities along the
%    boundary. The solution method is finite differencing on
%    a staggered grid with implicit diffusion and a Chorin
%    projection method for the pressure.
%    Visualization is done by a colormap-isoline plot for
%    pressure and normalized quiver and streamline plot for
%    the velocity field.
%    The standard setup solves a lid driven cavity problem.

% 07/2007 by Benjamin Seibold
%            http://www-math.mit.edu/~seibold/
% Feel free to modify for teaching and learning.
%-----------------------------------------------------------------------
switch example.case
    case 1 %Ë«¼ôÇÐ
        dt = 1e-4;    % time step
        tf = 1.8;    % final time
        lx = 1;       % width of box
        ly = 1;       % height of box
        nx = 128;      % number of x-gridpoints
        ny = 128;      % number of y-gridpoints 
        if example.subcase ==1
            %case 1
            delta = 0.05;
            rho = 1/30;
            Re = 1e4;     % Reynolds number
        else
        %case 2
            delta = 0.05;
            rho = 1/100;
            Re = 2e4;     % Reynolds number
        end
        % initial conditions
        x = linspace(0,lx,nx+1); hx = lx/nx;
        y = linspace(0,ly,ny+1); hy = ly/ny;
            x_0 = avg(x);
            y_0 = avg(y);
            x_1 = x(2:end-1);
            y_1 = y(2:end-1);
            idx = y_0<=0.5;
            U = (0*x_1)'+[tanh((y_0(idx)-0.25)/rho),tanh((0.75-y_0(~idx))/rho)];
            V = delta*sin(2*pi*x_0)'+0*y_1;
            % boundary conditions
            uW = [tanh((y_0(idx)-0.25)/rho),tanh((0.75-y_0(~idx))/rho)];
            uE = uW; 
            uN = 0*x+tanh(-0.25/rho);
            uS = 0*x+tanh(-0.25/rho);
            vN =  delta*sin(2*pi*x_0);
             vS = vN;
             vW = delta*sin(2*pi*0)+0*y;
            vE =  delta*sin(2*pi*1)+0*y;
    case  2 % Back_step
        %-----------------------------------------------------------------------
        Re = 50;     % Reynolds number
        dt = 0.01;    % time step
        tf = 20;    % final time
        lx = 20;       % width of box
        ly = 2;       % height of box
        nx = 60;      % number of x-gridpoints
        ny = 40;      % number of y-gridpoints

        IC = 1;
        x = linspace(0,lx,nx+1); hx = lx/nx;
        y = linspace(0,ly,ny+1); hy = ly/ny;
        %-----------------------------------------------------------------------
        % initial conditions
        % construct initial conditions
        if example.subcase ==1
              % parabolic initial conditions
            [u0,v0] = IC_parabolic(nx,ny);
        else
             % parabolic initial conditions with interpolation
            [u0,v0] = IC_parabolic_iterpolate(nx,ny);
        end
      
 
        
        U = u0(2:end-1,2:end-1)';
        V = v0(2:end-1,2:end-1)';
        % boundary conditions

        uN = u0(end,:);
        vN = v0(end,2:end-1);
        uS = u0(1,:);
        vS = v0(1,2:end-1);
        uW = u0(2:end-1,1)';
        vW = v0(:,1)';
        uE = U(end,:);
        vE = [0,V(end,:),0];

  
end


nsteps = 10;  % number of steps with graphic output
nt = ceil(tf/dt); dt = tf/nt;

%-----------------------------------------------------------------------
Ubc = dt/Re*([2*uS(2:end-1)' zeros(nx-1,ny-2) 2*uN(2:end-1)']/hy^2+...
      [uW;zeros(nx-3,ny);uE]/hx^2);
Vbc = dt/Re*([vS' zeros(nx,ny-3) vN']/hy^2+...
      [2*vW(2:end-1);zeros(nx-2,ny-1);2*vE(2:end-1)]/hx^2);






fprintf('initialization')
Lp = kron(speye(ny),K1(nx,hx,1))+kron(K1(ny,hy,1),speye(nx));
Lp(1,1) = 3/2*Lp(1,1);
perp = symamd(Lp); Rp = chol(Lp(perp,perp)); Rpt = Rp';
Lu = speye((nx-1)*ny)+dt/Re*(kron(speye(ny),K1(nx-1,hx,2))+...
     kron(K1(ny,hy,3),speye(nx-1)));
peru = symamd(Lu); Ru = chol(Lu(peru,peru)); Rut = Ru';
Lv = speye(nx*(ny-1))+dt/Re*(kron(speye(ny-1),K1(nx,hx,3))+...
     kron(K1(ny-1,hy,2),speye(nx)));
perv = symamd(Lv); Rv = chol(Lv(perv,perv)); Rvt = Rv';


fprintf(', time loop\n--20%%--40%%--60%%--80%%-100%%\n')
for k = 1:nt
   % treat nonlinear terms
   gamma = min(1.2*dt*max(max(max(abs(U)))/hx,max(max(abs(V)))/hy),1);
   Ue = [uW;U;uE]; Ue = [2*uS'-Ue(:,1) Ue 2*uN'-Ue(:,end)];
   Ve = [vS' V vN']; Ve = [2*vW-Ve(1,:);Ve;2*vE-Ve(end,:)];
   Ua = avg(Ue')'; Ud = diff(Ue')'/2;
   Va = avg(Ve);   Vd = diff(Ve)/2;
   UVx = diff(Ua.*Va-gamma*abs(Ua).*Vd)/hx;
   UVy = diff((Ua.*Va-gamma*Ud.*abs(Va))')'/hy;
   Ua = avg(Ue(:,2:end-1));   Ud = diff(Ue(:,2:end-1))/2;
   Va = avg(Ve(2:end-1,:)')'; Vd = diff(Ve(2:end-1,:)')'/2;
   U2x = diff(Ua.^2-gamma*abs(Ua).*Ud)/hx;
   V2y = diff((Va.^2-gamma*abs(Va).*Vd)')'/hy;
   U = U-dt*(UVy(2:end-1,:)+U2x);
   V = V-dt*(UVx(:,2:end-1)+V2y);
   
   % implicit viscosity
   rhs = reshape(U+Ubc,[],1);
   u(peru) = Ru\(Rut\rhs(peru));
   U = reshape(u,nx-1,ny);
   rhs = reshape(V+Vbc,[],1);
   v(perv) = Rv\(Rvt\rhs(perv));
   V = reshape(v,nx,ny-1);
   
   % pressure correction
   rhs = reshape(diff([uW;U;uE])/hx+diff([vS' V vN']')'/hy,[],1);
   p(perp) = -Rp\(Rpt\rhs(perp));
   P = reshape(p,nx,ny);
   U = U-diff(P)/hx;
   V = V-diff(P')'/hy;
   
   % boundary conditions
   if example.case == 1
    uW = 0.5*(U(1,:)+U(end,:));
    uE = uW; 
    uN = 0.5*[uW(1)+uW(end);(U(:,1)+U(:,end));uW(1)+uW(end)]';
    uS = uN;
    vN = 0.5*(V(:,1)+V(:,end))';
    vS = vN;
    vW =  0.5*[vN(1)+vN(end),(V(1,:)+V(end,:)),vN(1)+vN(end)];
    vE = vW;
     %-----------------------------------------------------------------------
    Ubc = dt/Re*([2*uS(2:end-1)' zeros(nx-1,ny-2) 2*uN(2:end-1)']/hx^2+...
          [uW;zeros(nx-3,ny);uE]/hy^2);
    Vbc = dt/Re*([vS' zeros(nx,ny-3) vN']/hx^2+...
          [2*vW(2:end-1);zeros(nx-2,ny-1);2*vE(2:end-1)]/hy^2);
   else
    uE = U(end,:); 
    vE = [0,V(end,:),0];
   end
   

   % visualization
   if floor(25*k/nt)>floor(25*(k-1)/nt), fprintf('.'), end
   if k==1|floor(nsteps*k/nt)>floor(nsteps*(k-1)/nt)
       if example.case == 1
          % stream function
          W = diff(U')'/hy-diff(V)/hx;
         % contourf(x(2:end-1),y(2:end-1),W',20,'w-'), hold on
          contour(x(2:end-1),y(2:end-1),W',25)
          title(sprintf('dt =  %0.2g nx = %d ny = %d Re = %0.1g   T = %0.2g',dt, nx, ny ,Re,k*dt))
       else
            Ue = [uS' avg([uW;U;uE]')' uN']';
            Ve = [vW;avg([vS' V vN']);vE]';   
            clf
            subplot(2,1,1)
            [X,Y] = meshgrid(0:hx:lx,0:hy:ly);
            quiver(X,Y,Ue,Ve)
            hold off,axis equal,axis([0 lx 0 ly])
            title(['vector plot, Re = ', num2str(Re)])       
            subplot(2,1,2)
            starty = 0:hy:ly;
            startx = starty*0;
            streamline(X,Y,Ue,Ve,startx,starty)
            hold off, axis equal,axis([0 lx 0 ly])
            suptitle(sprintf('hx = %0.2f hy = %0.2f Re = %0.1g   t = %1.3g',hx, hy, Re,k*dt))
       end
            drawnow
   end
    
    
end
fprintf('\n')

%=======================================================================

function B = avg(A,k)
if nargin<2, k = 1; end
if size(A,1)==1, A = A'; end
if k<2, B = (A(2:end,:)+A(1:end-1,:))/2; else, B = avg(A,k-1); end
if size(A,2)==1, B = B'; end

function A = K1(n,h,a11)
% a11: Neumann=1, Dirichlet=2, Dirichlet mid=3;
     A = spdiags([-1 a11 0;ones(n-2,1)*[-1 2 -1];0 a11 -1],-1:1,n,n)'/h^2;
