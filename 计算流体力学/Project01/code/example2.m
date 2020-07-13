function example = example2
%% This is the pde data with propagation of  Highly discontinuous data p71

T = 8;
IT = [-1,1];
meshize= struct('h',2/500,'tau',1.6/500,'IT',IT,'T',T);
option = {'upwind';'LF';'LW'};
pde = struct('exactu',@exactu,'uo',@uo,'a',1,'bdtype','Periodic','lp',2,'ordershow','ordershow','figureshow','figureshow');

example = struct('mesh',meshize,'option',{option},'pde',pde);
    function s = uo(x)
      L = length(x);
      s = zeros(1,L);
      xi = zeros(1,L);
      idxL = x<-0.7;
      idxR = (x>=-0.7);
      xi(idxL) = x(idxL)-0.3+2;
      xi(idxR) = x(idxR) - 0.3;
      idx1 = (xi<=-1/3);
      idx2 = logical((xi>-1/3).*(xi<1/3));
      idx3 = (xi>=1/3);
      xi1 = xi(idx1);
      xi2 = xi(idx2);
      xi3 = xi(idx3);
      s(idx1) = -xi1.*sin(1.5*pi*xi1.^2);
      s(idx2) = abs(sin(2*pi*xi2));
      s(idx3) = 2*xi3-1-sin(3*pi*xi3)/6;
    end

    function s = exactu(x,t)
       xo = mod(x - t+1,2)-1;
       s = uo(xo);
    end
        
        
end