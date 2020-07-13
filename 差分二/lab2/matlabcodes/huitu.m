R = 2;
theta =(0:-pi/3:-5*pi/3)';


x = R*[sin(theta),cos(theta)];
o = 0.3*(x(1:end,:)+[x(2:end,:);x(1,:)]);
oc = (x(1:end,:)+[x(2:end,:);x(1,:)]);