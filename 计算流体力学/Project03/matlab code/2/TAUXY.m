function tau_xy = TAUXY(mu,u,v,Dx,Dy,state,var)

switch state  % forecast or correct
    case 'forecast'
            switch var %E3 or E5 or F2 or F5
                case 'E3andE5'
                    u_y = Dyscheme(u,Dy,'center');
                    v_x =  Dxscheme(v,Dx,'backward');




                case 'F2andF5'
                    u_y = Dyscheme(u,Dy,'backward');
                    v_x =  Dxscheme(v,Dx,'center');




            end
    case 'correct'
            switch var
                case 'E3andE5'
                    u_y = Dyscheme(u,Dy,'center');
                    v_x =  Dxscheme(v,Dx,'forward');




                case 'F2andF5'
                    u_y = Dyscheme(u,Dy,'forward');
                    v_x =  Dxscheme(v,Dx,'center');



            end
        
end
        tau_xy = mu.*(u_y+v_x);
end








        
