function s = Partial_xmax_delta(x,y)
%% the partial_x of max_delta function
global delta;
s = 0.5*(1+(x-y)./sqrt((x-y).^2+delta^2));
end