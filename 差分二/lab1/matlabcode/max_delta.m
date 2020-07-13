function S = max_delta(x,y)
global delta;
S = 0.5*(x+y+sqrt((x-y).^2+delta^2));
end