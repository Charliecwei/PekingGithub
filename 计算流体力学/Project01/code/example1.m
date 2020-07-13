function example = example1
%% This is the pde data with propagation of a wave packet p69

T = 8;
IT = [T,T+1];
%IT = [0,1];
meshize= struct('h',1/200,'tau',0.8/200,'IT',IT,'T',T);
option = {'upwind';'LF';'LW'};
%option = {'LW'};
pde = struct('exactu',@exactu,'uo',@uo,'a',1,'bdtype','Fixed','lp',2,'ordershow','ordershow','figureshow','figureshow');

example = struct('mesh',meshize,'option',{option},'pde',pde);
    function s = uo(x)
           s= (exp(-100*(x-0.5).^2).*sin(80*x));
    end

    function s = exactu(x,t)
        s = uo(x-t);
    end
        
        
end