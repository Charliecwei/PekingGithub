function Is = CONVER(ep)
 %Check convergence
 err = max(abs(ep(:)))
 if err <=1e-8
     Is = 1;
 else
     Is = 0;
 end
end