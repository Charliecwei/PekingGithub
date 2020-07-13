function example = example4
%% This is the pde data with propagation of  Highly discontinuous data p71

T = 1.2;
IT = [-100,100];
Iu = [-3,3]; %A priori estimate of u
A = max(abs(fx(Iu(1):(Iu(2)-Iu(1))/100:Iu(2))));
Io = [IT(1)-A*T-5,IT(2)+A*T+5];

h = (Io(2)-Io(1))/1000;
tau = 0.25*h;
meshize= struct('h',h,'tau',tau,'IT',IT,'T',T);
option = {'upwind';'LF';'LW'};
a = struct('f',@f,'a',@fx,'Iu',Iu);
pde = struct('uo',@uo,'a',a,'bdtype','Fixed');








example = struct('mesh',meshize,'option',{option},'pde',pde);
    function s = uo(x)
      s =  -2*sign(x);
    end

 
function s = f(u)
   s = 0.25*(u.^2-1).*(u.^2-4);
end

function  s = fx(u)
 s= u.^3-2.5*u;
end
        
        
end