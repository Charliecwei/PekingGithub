function example = example3
%% This is the pde data with propagation of  Highly discontinuous data p71

T = 2;
IT = [0,2*pi];
h = (IT(2)-IT(1))/100;
tau = 0.5*h;
meshize= struct('h',h,'tau',tau,'IT',IT,'T',T);
option = {'upwind';'LF';'LW'};
a = struct('f',@fburger,'a',@fxburger);
pde = struct('uo',@uo,'a',a,'bdtype','Periodic');


example = struct('mesh',meshize,'option',{option},'pde',pde);
    function s = uo(x)
      s =  0.5+sin(x);
    end

 
function s = fburger(u)
   s = 0.5*u.^2;
end

function  s = fxburger(u)
 s= u;
end
        
        
end