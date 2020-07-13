function tau_xx = TAUXX(mu,u,v,Dx,Dy,state)
lambda = -2*mu/3;
switch state
    case 'forecast'
                u_x  = Dxscheme(u,Dx,'backward');
                v_y  = Dyscheme(v,Dy,'center');
                
    case 'correct'
               u_x  = Dxscheme(u,Dx,'forward');
               v_y  = Dyscheme(v,Dy,'center');
               
end

tau_xx = lambda.*(u_x+v_y)+2*mu.*u_x;
end
        
