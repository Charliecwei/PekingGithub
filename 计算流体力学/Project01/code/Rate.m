function [rate,err,hs] = Rate(pde,var,T,h,tau,Io,bdtype,lp)
%% compute the rate fo var method at time T,tau/h is pde.a constant for lp norm 1<=lp<inf
x = Io(1):h:Io(2);
lam = tau/h;
rate = zeros(1,4);
err = zeros(1,5);
hs = zeros(1,5);
hs(1) = h;
Uo = pde.exactu(x,0);
  switch var
      case 'upwind'
          [U, CourantNum]= UpWindCIR(Uo,h,tau,pde.a,T,bdtype);
      case 'LF'
          [U,CourantNum] = LaxFriedrichs(Uo,h,tau,pde.a,T,bdtype);
      case 'LW'
          [U,CourantNum] = LaxWendroff(Uo,h,tau,pde.a,T,bdtype);
  end

  switch bdtype %compute err
      case 'Periodic'
          if isinf(lp)
            err(1) = norm(U-pde.exactu(x,T),lp);
          else
              err(1) = power(h,1/lp)*norm(U-pde.exactu(x,T),lp);
          end
      case 'Fixed'
              as = CourantNum/lam;
              IT=[Io(1)+as*T+1,Io(2)-as*T-1];
              idx = (x<=IT(2))&(x>=IT(1));
              if isinf(lp)
                     err(1) = norm(U(idx)-pde.exactu(x(idx),T),lp);
              else
                    err(1) = power(h,1/lp)*norm(U(idx)-pde.exactu(x(idx),T),lp);
              end
  end
 
  
  for i = 1:4
      h = h/2;
      hs(i+1) = h;
      x = Io(1):h:Io(2);
      tau = lam*h;
      Uo = pde.exactu(x,0);    
   switch var
      case 'upwind'
          [U, CourantNum]= UpWindCIR(Uo,h,tau,pde.a,T,bdtype);
      case 'LF'
          [U,CourantNum] = LaxFriedrichs(Uo,h,tau,pde.a,T,bdtype);
      case 'LW'
          [U,CourantNum] = LaxWendroff(Uo,h,tau,pde.a,T,bdtype);
   end
  
    switch bdtype %compute err
      case 'Periodic'
          if isinf(lp)
            err(i+1) = norm(U-pde.exactu(x,T),lp);
          else
              err(i+1) = power(h,1/lp)*norm(U-pde.exactu(x,T),lp);
          end
      case 'Fixed'
              as = CourantNum/lam;
              IT=[Io(1)+as*T+1,Io(2)-as*T-1];
              idx = (x<=IT(2))&(x>=IT(1));
              if isinf(lp)
                     err(i+1) = norm(U(idx)-pde.exactu(x(idx),T),lp);
              else
                    err(i+1) = power(h,1/lp)*norm(U(idx)-pde.exactu(x(idx),T),lp);
              end
    end
  end
  
  
  %% compute rate
  for i = 1:4
      rate(i) = (log(err(i))-log(err(i+1)))/log(2);
  end
      
            