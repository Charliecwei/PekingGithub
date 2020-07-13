function s = Partial_xmin_delta(x,y)
%% the partial_x of min_delta function
global  delta;
s = 0.5*(1-(x-y)./sqrt((x-y).^2+delta^2));
end