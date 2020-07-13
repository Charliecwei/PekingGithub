function tau_yy = TAUYY(mu,u,v,Dx,Dy,state)
lambda = -2*mu/3;
switch state
    case 'forecast'
                u_x  = Dxscheme(u,Dx,'center');
                v_y  = Dyscheme(v,Dy,'backward');
                
    case 'correct'
               u_x  = Dxscheme(u,Dx,'center');
               v_y  = Dyscheme(v,Dy,'forward');
               
end
%tau_yy = v_y;
tau_yy = lambda.*(u_x+v_y)+2*mu.*v_y;
end
     